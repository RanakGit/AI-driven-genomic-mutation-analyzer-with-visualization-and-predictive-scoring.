import csv
import gzip
import os
import re
from urllib.parse import urlparse


def read_clinvar_tsv(file_input):
    """Read simple ClinVar tab-delimited variant file and return normalized list."""
    if hasattr(file_input, 'read'):
        raw = file_input.read()
        if isinstance(raw, bytes):
            text = raw.decode('utf-8', errors='ignore')
        else:
            text = raw
        try:
            file_input.seek(0)
        except Exception:
            pass
        lines = text.splitlines()
    elif isinstance(file_input, str):
        with open(file_input, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.read().splitlines()
    else:
        raise TypeError('Input must be file path or file-like object')

    reader = csv.DictReader(lines, delimiter='\t')
    variants = []
    for row in reader:
        variants.append({
            'chrom': row.get('Chromosome', row.get('CHROM', '')).strip(),
            'pos': int(row.get('Pos', row.get('POS', 0))) if row.get('Pos', row.get('POS', '')).strip().isdigit() else None,
            'ref': row.get('Ref', row.get('REF', '')).strip(),
            'alt': row.get('Alt', row.get('ALT', '')).strip(),
            'clinical_significance': row.get('ClinicalSignificance', row.get('CLNSIG', '')).strip()
        })
    return variants


def download_clinvar_vcf(url, target_path):
    """Download ClinVar VCF file from URL and save to target_path."""
    try:
        import requests
    except ImportError as e:
        raise RuntimeError('requests is required for downloading ClinVar data: pip install requests') from e

    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()

    if urlparse(url).path.endswith('.gz') or target_path.endswith('.gz'):
        with open(target_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
    else:
        with open(target_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

    return target_path


def parse_clinvar_vcf(file_input):
    """Parse ClinVar VCF and return list of normalized variant dictionaries."""
    def open_file(path_or_obj):
        if hasattr(path_or_obj, 'read'):
            return path_or_obj
        if str(path_or_obj).endswith('.gz'):
            return gzip.open(path_or_obj, 'rt', encoding='utf-8', errors='ignore')
        return open(path_or_obj, 'r', encoding='utf-8', errors='ignore')

    reader = open_file(file_input)
    variants = []
    for line in reader:
        if not line or line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue
        chrom, pos, id_, ref, alt, qual, filt, info = parts[:8]
        cs = 'Unknown'
        match = re.search(r'CLNSIG=([^;]+)', info)
        if match:
            cs = match.group(1).split('|')[0]

        variants.append({
            'chrom': chrom,
            'pos': int(pos),
            'ref': ref,
            'alt': alt,
            'clinical_significance': cs,
            'id': id_
        })

    if not hasattr(file_input, 'read'):
        if not str(file_input).endswith('.gz') and hasattr(reader, 'close'):
            reader.close()
    return variants


def variants_to_training_records(vcf_variants, ref_seq=''):
    """Convert raw VCF variants to training records with features for ML."""
    records = []
    for v in vcf_variants:
        vtype = 'SNP'
        if len(v['ref']) > 1 or len(v['alt']) > 1:
            if len(v['ref']) == 1 and len(v['alt']) > 1:
                vtype = 'INS'
            elif len(v['alt']) == 1 and len(v['ref']) > 1:
                vtype = 'DEL'
            else:
                vtype = 'MNV'

        aa_change = f"{v['ref']}>{v['alt']}"
        cs = v.get('clinical_significance', 'Unknown')

        record = {
            'type': vtype,
            'aa_change': aa_change,
            'clinical_significance': cs,
            'pos': v.get('pos'),
            'chrom': v.get('chrom'),
            'id': v.get('id')
        }

        if ref_seq and v.get('pos') is not None and vtype == 'SNP':
            pc = (v['pos'] - 1) // 3
            record['codon_index'] = pc
            codon = ref_seq[pc*3:(pc*3)+3]
            record['codon_ref'] = codon
            if len(codon) == 3:
                record['codon_change'] = f"{codon}:{v['ref']}>{v['alt']}"

        records.append(record)
    return records
