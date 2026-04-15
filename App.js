import React, { useState } from 'react';
import axios from 'axios';

function App() {
  const [referenceSeq, setReferenceSeq] = useState('');
  const [querySeq, setQuerySeq] = useState('');
  const [result, setResult] = useState(null);
  const [error, setError] = useState('');

  const runAnalysis = async () => {
    setError('');
    try {
      const resp = await axios.post('/api/analyze', { reference_seq: referenceSeq, query_seq: querySeq });
      setResult(resp.data.result);
    } catch (e) {
      setError(e.response?.data?.message || e.message);
    }
  };

  return (
    <div style={{ padding: 24, fontFamily: 'Arial, sans-serif' }}>
      <h1>SmartVariant Frontend</h1>
      <div>
        <label>Reference sequence:</label><br />
        <textarea value={referenceSeq} onChange={(e) => setReferenceSeq(e.target.value)} rows={3} cols={80}></textarea>
      </div>
      <div>
        <label>Query sequence:</label><br />
        <textarea value={querySeq} onChange={(e) => setQuerySeq(e.target.value)} rows={3} cols={80}></textarea>
      </div>
      <button onClick={runAnalysis}>Analyze variants</button>
      {error && <div style={{ color: 'red' }}>{error}</div>}
      {result && (
        <div style={{ marginTop: 20 }}>
          <h2>Result</h2>
          <p>Risk score: {result.risk_score}</p>
          <pre>{JSON.stringify(result, null, 2)}</pre>
        </div>
      )}
    </div>
  );
}

export default App;
