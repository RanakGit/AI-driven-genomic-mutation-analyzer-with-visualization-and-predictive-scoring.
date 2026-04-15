import os
from flask import Flask, request, redirect, url_for, render_template, flash, jsonify
import plotly
import plotly.graph_objs as go
import json

from utils.parser import read_sequence
from utils.mutation import detect_variants, summarize_variants
from utils.database import init_db, save_analysis, get_history
from utils.clinvar import read_clinvar_tsv, download_clinvar_vcf, parse_clinvar_vcf, variants_to_training_records
from utils.protein_effect import variant_protein_impact
from model.impact_predictor import load_or_train_model, train_clinvar_model, predict_variant_impact, load_model_info


app = Flask(__name__)
app.secret_key = 'smartvariant_secret_key'

# Initialize DB and ML model once
init_db()
pipeline = load_or_train_model()


def analyze_pair(ref_seq, query_seq):
    base_analysis = detect_variants(ref_seq, query_seq)
    variants = base_analysis['variants']
    counts = summarize_variants(variants)

    # Predict variant impacts
    impacts = []
    impact_scores = {
        'Benign': 0,
        'Harmful': 1
    }
    tot_score = 0
    for v in variants:
        ml_impact = predict_variant_impact(pipeline, v)
        protein_impact = variant_protein_impact(ref_seq, v)
        combined_impact = ml_impact if ml_impact != 'Unknown' else protein_impact
        impacts.append({
            **v,
            'ml_impact': ml_impact,
            'protein_impact': protein_impact,
            'impact': combined_impact
        })
        tot_score += impact_scores.get(ml_impact, 0)

    risk_score = float(tot_score / max(1, len(variants)))
    return {
        'variants': impacts,
        'summary': counts,
        'risk_score': risk_score,
        'alignment_score': base_analysis['score']
    }


@app.route('/', methods=['GET'])
def index():
    history = get_history(limit=20)
    return render_template('index.html', history=history)


@app.route('/analyze', methods=['POST'])
def analyze():
    if 'reference_seq' not in request.files or 'query_seq' not in request.files:
        flash('Both reference and query sequence files are required')
        return redirect(url_for('index'))

    ref_file = request.files['reference_seq']
    query_file = request.files['query_seq']

    # progress indicator can be managed in JS only for asynchronous calls.
    # For now, template-level messages can show status updates.
    flash('Analysis is running... this may take several seconds.')

    if ref_file.filename == '' or query_file.filename == '':
        flash('Both sequence files must be selected.')
        return redirect(url_for('index'))

    try:
        ref_records = read_sequence(ref_file)
        query_records = read_sequence(query_file)

        if not ref_records or not query_records:
            raise ValueError('No sequences found in uploaded files')

        ref_entry = ref_records[0]
        query_entry = query_records[0]

        result = analyze_pair(ref_entry['sequence'], query_entry['sequence'])

        save_analysis(
            ref_entry.get('id', 'ref'),
            query_entry.get('id', 'query'),
            result['summary'].get('SNP', 0),
            result['summary'].get('INS', 0),
            result['summary'].get('DEL', 0),
            result['risk_score']
        )

        # build plotly bar chart for mutation distribution
        bar = go.Bar(
            x=['SNP', 'INS', 'DEL'],
            y=[result['summary'].get('SNP', 0), result['summary'].get('INS', 0), result['summary'].get('DEL', 0)],
            marker=dict(color=['blue', 'orange', 'red'])
        )
        chart = json.dumps([bar], cls=plotly.utils.PlotlyJSONEncoder)

        return render_template('results.html',
                               reference_id=ref_entry.get('id', 'ref'),
                               query_id=query_entry.get('id', 'query'),
                               result=result,
                               chart=chart)

    except Exception as ex:
        flash(str(ex))
        return redirect(url_for('index'))


@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    payload = request.get_json(silent=True)
    if payload and 'reference_seq' in payload and 'query_seq' in payload:
        ref_seq = payload['reference_seq']
        query_seq = payload['query_seq']
    else:
        return jsonify({'error': 'Must provide JSON with reference_seq and query_seq'}), 400

    try:
        result = analyze_pair(ref_seq, query_seq)
        return jsonify({'status': 'success', 'result': result})
    except Exception as ex:
        return jsonify({'status': 'error', 'message': str(ex)}), 500


@app.route('/api/batch', methods=['POST'])
def api_batch():
    payload = request.get_json(silent=True)
    if not payload or 'pairs' not in payload or not isinstance(payload['pairs'], list):
        return jsonify({'error': 'Must provide JSON with "pairs" list of objects.'}), 400

    outputs = []
    for item in payload['pairs']:
        ref_seq = item.get('reference_seq')
        query_seq = item.get('query_seq')
        if not ref_seq or not query_seq:
            outputs.append({'error': 'Missing reference_seq or query_seq in pair'})
            continue
        try:
            analysis = analyze_pair(ref_seq, query_seq)
            outputs.append({'status': 'success', 'result': analysis})
        except Exception as ex:
            outputs.append({'status': 'error', 'message': str(ex)})

    return jsonify({'batch_results': outputs})


@app.route('/api/clinvar', methods=['POST'])
def api_clinvar():
    if 'clinvar_file' not in request.files:
        return jsonify({'error': 'No clinvar_file uploaded'}), 400

    file = request.files['clinvar_file']
    try:
        variants = read_clinvar_tsv(file)
        global pipeline
        pipeline = train_clinvar_model(variants)

        model_info = load_model_info()
        return jsonify({'status': 'success', 'source': 'tsv', 'count': len(variants), 'model_info': model_info})
    except Exception as ex:
        return jsonify({'status': 'error', 'message': str(ex)}), 500


@app.route('/api/clinvar/download', methods=['POST'])
def api_clinvar_download():
    payload = request.get_json(silent=True) or {}
    url = payload.get('url')
    algorithm = payload.get('algorithm', 'rf')

    if not url:
        return jsonify({'status': 'error', 'message': 'URL is required for ClinVar download'}), 400

    try:
        tmp_path = os.path.join('data', 'clinvar_download.vcf.gz')
        download_clinvar_vcf(url, tmp_path)
        vcf_variants = parse_clinvar_vcf(tmp_path)
        training_records = variants_to_training_records(vcf_variants)
        global pipeline
        pipeline = train_clinvar_model(training_records, algorithm=algorithm)
        model_info = load_model_info()
        return jsonify({'status': 'success', 'downloaded': len(vcf_variants), 'trained': len(training_records), 'model_info': model_info})
    except Exception as ex:
        return jsonify({'status': 'error', 'message': str(ex)}), 500


@app.route('/history')
def history():
    rows = get_history(limit=100)
    return render_template('history.html', history=rows)


@app.route('/metrics')
def metrics():
    rows = get_history(limit=200)
    # metric: per run risk over time + mutation load
    runs = []
    for r in reversed(rows):
        runs.append({
            'id': r[0],
            'timestamp': r[1],
            'snp': r[4],
            'ins': r[5],
            'del': r[6],
            'risk': float(r[7] or 0)
        })

    return render_template('metrics.html', runs=runs)


if __name__ == '__main__':
    app.run(debug=True)

