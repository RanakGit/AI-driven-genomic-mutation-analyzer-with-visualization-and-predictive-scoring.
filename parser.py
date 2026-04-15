import io
from Bio import SeqIO


def read_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append({
            "id": record.id,
            "sequence": str(record.seq),
            "length": len(record.seq)
        })
    return sequences


def read_genbank(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "genbank"):
        sequences.append({
            "id": record.id,
            "name": record.name,
            "description": record.description,
            "sequence": str(record.seq),
            "length": len(record.seq)
        })
    return sequences


def _is_path(obj):
    return isinstance(obj, str)


def _ensure_text_stream(stream):
    if hasattr(stream, 'read'):
        try:
            sample = stream.read(0)
            if isinstance(sample, bytes):
                stream = io.TextIOWrapper(stream, encoding='utf-8', newline='')
        except Exception:
            pass
    return stream


def _build_records(records):
    return [
        {
            'id': rec.id,
            'sequence': str(rec.seq),
            'length': len(rec.seq)
        }
        for rec in records
    ]


def read_sequence(file_input):
    if _is_path(file_input):
        file_name = file_input
        lower_name = file_name.lower()
        if lower_name.endswith(('.fasta', '.fa', '.fna')):
            return read_fasta(file_name)
        elif lower_name.endswith(('.gb', '.gbk')):
            return read_genbank(file_name)
        else:
            raise ValueError(f"Unsupported file extension: {file_name}")

    filename = getattr(file_input, 'filename', getattr(file_input, 'name', ''))
    lower_name = filename.lower() if filename else ''
    stream = getattr(file_input, 'stream', file_input)
    stream = _ensure_text_stream(stream)
    stream.seek(0)

    if lower_name.endswith(('.fasta', '.fa', '.fna')):
        fmt = 'fasta'
    elif lower_name.endswith(('.gb', '.gbk')):
        fmt = 'genbank'
    else:
        records = list(SeqIO.parse(stream, 'fasta'))
        if records:
            stream.seek(0)
            return _build_records(SeqIO.parse(stream, 'fasta'))
        stream.seek(0)
        return _build_records(SeqIO.parse(stream, 'genbank'))

    stream.seek(0)
    return _build_records(SeqIO.parse(stream, fmt))

