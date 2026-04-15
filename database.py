import sqlite3
from datetime import datetime

DB_PATH = 'data/analysis_history.db'


def init_db():
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS analysis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp TEXT NOT NULL,
            ref_id TEXT,
            query_id TEXT,
            snp_count INTEGER,
            ins_count INTEGER,
            del_count INTEGER,
            risk_score REAL
        )
    ''')
    conn.commit()
    conn.close()


def save_analysis(ref_id, query_id, snp_count, ins_count, del_count, risk_score):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute('''
        INSERT INTO analysis (timestamp, ref_id, query_id, snp_count, ins_count, del_count, risk_score)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    ''', (datetime.utcnow().isoformat(), ref_id, query_id, snp_count, ins_count, del_count, risk_score))
    conn.commit()
    conn.close()


def get_history(limit=100):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute('''
        SELECT id, timestamp, ref_id, query_id, snp_count, ins_count, del_count, risk_score
        FROM analysis
        ORDER BY id DESC
        LIMIT ?
    ''', (limit,))
    rows = c.fetchall()
    conn.close()
    return rows
