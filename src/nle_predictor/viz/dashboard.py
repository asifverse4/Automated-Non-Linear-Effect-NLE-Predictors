"""
Interactive Standalone HTML Dashboard Generator for NLE-Predictor.
Produces a self-contained, offline-compatible HTML application with
real-time Newton-Raphson solvers, dynamic charts, preset benchmarks,
and CSV overlay capabilities.
"""

import os


def generate_interactive_html(output_path: str = "manuscript/interactive_dashboard.html") -> str:
    """Generate self-contained interactive HTML dashboard."""
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>NLE-Predictor: Interactive Chemistry Dashboard</title>
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>
  <style>
    :root {
      --primary: #1e3a8a;
      --primary-light: #3b82f6;
      --secondary: #0d9488;
      --accent: #dc2626;
      --bg: #f8fafc;
      --card-bg: #ffffff;
      --text: #0f172a;
      --text-muted: #64748b;
      --border: #e2e8f0;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; }
    body { background: var(--bg); color: var(--text); padding: 24px; line-height: 1.5; }
    .header { text-align: center; margin-bottom: 24px; }
    .header h1 { font-size: 1.85rem; color: var(--primary); margin-bottom: 6px; }
    .header p { color: var(--text-muted); font-size: 0.95rem; }
    .container { max-width: 1300px; margin: 0 auto; display: grid; grid-template-columns: 360px 1fr; gap: 24px; }
    @media (max-width: 960px) { .container { grid-template-columns: 1fr; } }
    .card { background: var(--card-bg); border-radius: 12px; border: 1px solid var(--border); padding: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.05); }
    .card h2 { font-size: 1.15rem; color: var(--primary); margin-bottom: 16px; border-bottom: 2px solid var(--border); padding-bottom: 8px; }
    .control-group { margin-bottom: 16px; }
    .control-group label { display: flex; justify-content: space-between; font-weight: 500; font-size: 0.88rem; margin-bottom: 6px; }
    .control-group input[type="range"] { width: 100%; cursor: pointer; accent-color: var(--primary-light); }
    .presets { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 18px; }
    .preset-btn { background: #f1f5f9; border: 1px solid #cbd5e1; padding: 6px 10px; font-size: 0.8rem; border-radius: 6px; cursor: pointer; transition: all 0.2s; }
    .preset-btn:hover { background: var(--primary-light); color: white; border-color: var(--primary-light); }
    .charts-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 20px; }
    @media (max-width: 768px) { .charts-grid { grid-template-columns: 1fr; } }
    .chart-box { background: var(--card-bg); border-radius: 12px; border: 1px solid var(--border); padding: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.05); height: 380px; }
    .metrics { display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px; margin-top: 16px; }
    @media (max-width: 768px) { .metrics { grid-template-columns: repeat(2, 1fr); } }
    .metric-card { background: #f8fafc; border: 1px solid var(--border); border-radius: 8px; padding: 12px; text-align: center; }
    .metric-card .val { font-size: 1.35rem; font-weight: bold; color: var(--primary); margin-top: 4px; }
    .metric-card .label { font-size: 0.78rem; color: var(--text-muted); text-transform: uppercase; letter-spacing: 0.5px; }
    .upload-area { border: 2px dashed #cbd5e1; border-radius: 8px; padding: 12px; text-align: center; cursor: pointer; margin-top: 12px; font-size: 0.85rem; color: var(--text-muted); }
    .upload-area:hover { border-color: var(--primary-light); color: var(--primary); }
    footer { text-align: center; margin-top: 32px; font-size: 0.85rem; color: var(--text-muted); }
  </style>
</head>
<body>
  <div class="header">
    <h1>NLE-Predictor: Interactive Chemistry Dashboard</h1>
    <p>Rigorous Thermodynamic and Kinetic Modeling Framework &bull; Department of Chemistry, Jamia Hamdard</p>
  </div>

  <div class="container">
    <!-- Controls Column -->
    <div class="card">
      <h2>Model Controls & Parameters</h2>

      <div class="control-group">
        <label>Benchmark Presets</label>
        <div class="presets">
          <button class="preset-btn" onclick="applyPreset(-1.74, 1e5, 1.0, 0.0)">Kagan 1986</button>
          <button class="preset-btn" onclick="applyPreset(-3.5, 1e6, 0.95, 0.0)">Noyori DAIB</button>
          <button class="preset-btn" onclick="applyPreset(-1.5, 5e4, 0.98, 0.0)">Walsh BINOL</button>
          <button class="preset-btn" onclick="applyPreset(-0.41, 1e5, 1.0, 0.0)">Ideal Linear</button>
          <button class="preset-btn" onclick="applyPreset(0.38, 1e5, 1.0, 0.0)">(-)NLE Depletion</button>
        </div>
      </div>

      <div class="control-group">
        <label><span>&Delta;G (RS - RR):</span> <span id="val_dg">-1.80 kcal/mol</span></label>
        <input type="range" id="slider_dg" min="-4.0" max="2.0" step="0.05" value="-1.80">
      </div>

      <div class="control-group">
        <label><span>Association log10(K_homo):</span> <span id="val_khomo">5.00</span></label>
        <input type="range" id="slider_khomo" min="1.0" max="7.0" step="0.1" value="5.0">
      </div>

      <div class="control-group">
        <label><span>Intrinsic Enantioselectivity (ee0):</span> <span id="val_ee0">1.00</span></label>
        <input type="range" id="slider_ee0" min="0.1" max="1.0" step="0.01" value="1.00">
      </div>

      <div class="control-group">
        <label><span>Active Dimer Rate g (kD/kM):</span> <span id="val_g">0.00</span></label>
        <input type="range" id="slider_g" min="0.0" max="2.0" step="0.05" value="0.00">
      </div>

      <div class="control-group">
        <label><span>Temperature (T):</span> <span id="val_temp">298.15 K</span></label>
        <input type="range" id="slider_temp" min="200" max="350" step="1" value="298">
      </div>

      <div class="control-group">
        <label>Experimental Data Overlay (CSV)</label>
        <div class="upload-area" onclick="document.getElementById('csv_input').click()">
          Drag & Drop or Click to Upload CSV
          <input type="file" id="csv_input" accept=".csv" style="display: none;" onchange="handleCSV(event)">
        </div>
      </div>
    </div>

    <!-- Charts & Metrics Column -->
    <div>
      <div class="charts-grid">
        <div class="chart-box">
          <canvas id="chartNLE"></canvas>
        </div>
        <div class="chart-box">
          <canvas id="chartSpeciation"></canvas>
        </div>
      </div>

      <div class="card">
        <h2>Mathematical Curve Diagnostics</h2>
        <div class="metrics">
          <div class="metric-card">
            <div class="label">Classification</div>
            <div class="val" id="diag_class">(+)NLE</div>
          </div>
          <div class="metric-card">
            <div class="label">Initial Slope (S0)</div>
            <div class="val" id="diag_s0">8.90</div>
          </div>
          <div class="metric-card">
            <div class="label">Max Amplification</div>
            <div class="val" id="diag_max_amp">10.75</div>
          </div>
          <div class="metric-card">
            <div class="label">Area Under Asymmetry</div>
            <div class="val" id="diag_aua">0.379</div>
          </div>
        </div>
      </div>
    </div>
  </div>

  <footer>
    NLE-Predictor &bull; Open-Source Computational Instrument &bull; Author: Asif Raza (asifrazakne2005@gmail.com)
  </footer>

  <script>
    const R_GAS = 0.00198720425864083;
    let experimentalData = null;
    let chartNLEInstance = null;
    let chartSpecInstance = null;

    function solveMassBalance(cR, cS, kHomo, kHetero) {
      let mR = (-1 + Math.sqrt(Math.max(1e-12, 1 + 8 * kHomo * cR))) / (4 * kHomo);
      let mS = (-1 + Math.sqrt(Math.max(1e-12, 1 + 8 * kHomo * cS))) / (4 * kHomo);

      for (let iter = 0; iter < 40; iter++) {
        let f1 = mR + 2 * kHomo * mR * mR + kHetero * mR * mS - cR;
        let f2 = mS + 2 * kHomo * mS * mS + kHetero * mR * mS - cS;
        if (Math.abs(f1) < 1e-9 && Math.abs(f2) < 1e-9) break;

        let j11 = 1 + 4 * kHomo * mR + kHetero * mS;
        let j12 = kHetero * mR;
        let j21 = kHetero * mS;
        let j22 = 1 + 4 * kHomo * mS + kHetero * mR;
        let det = j11 * j22 - j12 * j21;
        if (Math.abs(det) < 1e-15) break;

        let dmR = (j22 * f1 - j12 * f2) / det;
        let dmS = (-j21 * f1 + j11 * f2) / det;
        mR = Math.max(0, mR - dmR);
        mS = Math.max(0, mS - dmS);
      }
      return [mR, mS];
    }

    function computeModel() {
      const dG = parseFloat(document.getElementById('slider_dg').value);
      const logK = parseFloat(document.getElementById('slider_khomo').value);
      const ee0 = parseFloat(document.getElementById('slider_ee0').value);
      const g = parseFloat(document.getElementById('slider_g').value);
      const temp = parseFloat(document.getElementById('slider_temp').value);

      document.getElementById('val_dg').innerText = dG.toFixed(2) + ' kcal/mol';
      document.getElementById('val_khomo').innerText = logK.toFixed(2);
      document.getElementById('val_ee0').innerText = ee0.toFixed(2);
      document.getElementById('val_g').innerText = g.toFixed(2);
      document.getElementById('val_temp').innerText = temp.toFixed(0) + ' K';

      const kHomo = Math.pow(10, logK);
      const kRatio = 2 * Math.exp(-dG / (R_GAS * temp));
      const kHetero = kRatio * kHomo;

      const nSteps = 101;
      const eeCatArr = [];
      const eeProdArr = [];
      const mRArr = [];
      const mSArr = [];
      const dRRArr = [];
      const dRSArr = [];

      for (let i = 0; i < nSteps; i++) {
        let eeCat = i / (nSteps - 1);
        eeCatArr.push(eeCat * 100);

        let cR = (1 + eeCat) / 2;
        let cS = (1 - eeCat) / 2;
        let [mR, mS] = solveMassBalance(cR, cS, kHomo, kHetero);

        let dRR = kHomo * mR * mR;
        let dRS = kHetero * mR * mS;

        mRArr.push(mR);
        mSArr.push(mS);
        dRRArr.push(dRR);
        dRSArr.push(dRS);

        let eeProd;
        if (g > 0) {
          let num = (mR - mS) * (ee0 + 2 * g * ee0 * kHomo * (mR + mS));
          let den = (mR + mS) + 2 * g * kHomo * (mR * mR + mS * mS) + 2 * g * kHetero * mR * mS;
          eeProd = den > 1e-12 ? num / den : ee0 * eeCat;
        } else {
          eeProd = (mR + mS) > 1e-12 ? ee0 * (mR - mS) / (mR + mS) : ee0 * eeCat;
        }
        eeProdArr.push(Math.max(-1, Math.min(1, eeProd)) * 100);
      }

      // Diagnostics
      let s0 = (eeProdArr[1] - eeProdArr[0]) / (eeCatArr[1] - eeCatArr[0]);
      let maxAmp = 1.0;
      for (let i = 1; i < nSteps - 1; i++) {
        let ratio = eeProdArr[i] / (eeCatArr[i] * ee0);
        if (ratio > maxAmp) maxAmp = ratio;
      }
      let aua = 0;
      for (let i = 0; i < nSteps - 1; i++) {
        let dx = (eeCatArr[i+1] - eeCatArr[i]) / 100;
        let yLinear = ee0 * ((eeCatArr[i+1] + eeCatArr[i]) / 2);
        let yModel = (eeProdArr[i+1] + eeProdArr[i]) / 2;
        aua += (yModel - yLinear) * dx;
      }

      document.getElementById('diag_s0').innerText = s0.toFixed(2);
      document.getElementById('diag_max_amp').innerText = maxAmp.toFixed(2);
      document.getElementById('diag_aua').innerText = (aua / 100).toFixed(4);
      let classification = "(+)NLE";
      if (Math.abs(dG - (-R_GAS * temp * Math.log(2))) < 0.1) classification = "Linear";
      else if (dG > -R_GAS * temp * Math.log(2)) classification = "(-)NLE";
      document.getElementById('diag_class').innerText = classification;

      updateCharts(eeCatArr, eeProdArr, mRArr, mSArr, dRRArr, dRSArr);
    }

    function updateCharts(eeCat, eeProd, mR, mS, dRR, dRS) {
      // Chart 1: NLE
      const datasetsNLE = [
        { label: 'NLE-Predictor Model', data: eeProd.map((y, i) => ({ x: eeCat[i], y: y })), borderColor: '#1e3a8a', backgroundColor: 'transparent', borderWidth: 2.5, pointRadius: 0, tension: 0.1 },
        { label: 'Ideal Linear Reference', data: eeCat.map(x => ({ x: x, y: x })), borderColor: '#94a3b8', borderDash: [5, 5], borderWidth: 1.5, pointRadius: 0 }
      ];

      if (experimentalData) {
        datasetsNLE.push({
          label: 'Experimental Data',
          data: experimentalData,
          backgroundColor: '#dc2626',
          borderColor: '#991b1b',
          pointRadius: 5,
          showLine: false,
          type: 'scatter'
        });
      }

      if (!chartNLEInstance) {
        const ctxNLE = document.getElementById('chartNLE').getContext('2d');
        chartNLEInstance = new Chart(ctxNLE, {
          type: 'line',
          data: { datasets: datasetsNLE },
          options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: { title: { display: true, text: 'Product Enantiopurity vs Catalyst Optical Purity' } },
            scales: {
              x: { type: 'linear', min: 0, max: 100, title: { display: true, text: 'Catalyst ee (%)' } },
              y: { min: 0, max: 100, title: { display: true, text: 'Product ee (%)' } }
            }
          }
        });
      } else {
        chartNLEInstance.data.datasets = datasetsNLE;
        chartNLEInstance.update('none');
      }

      // Chart 2: Speciation
      const datasetsSpec = [
        { label: '[MR] Major Monomer', data: mR.map((y, i) => ({ x: eeCat[i], y: y })), borderColor: '#2563eb', borderWidth: 2, pointRadius: 0 },
        { label: '[MS] Minor Monomer', data: mS.map((y, i) => ({ x: eeCat[i], y: y })), borderColor: '#9333ea', borderWidth: 2, pointRadius: 0 },
        { label: '[(R,S)] Heterochiral Dimer', data: dRS.map((y, i) => ({ x: eeCat[i], y: y })), borderColor: '#059669', borderWidth: 2, pointRadius: 0 },
        { label: '[(R,R)] Homochiral Dimer', data: dRR.map((y, i) => ({ x: eeCat[i], y: y })), borderColor: '#d97706', borderWidth: 2, pointRadius: 0 }
      ];

      if (!chartSpecInstance) {
        const ctxSpec = document.getElementById('chartSpeciation').getContext('2d');
        chartSpecInstance = new Chart(ctxSpec, {
          type: 'line',
          data: { datasets: datasetsSpec },
          options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: { title: { display: true, text: 'Equilibrium Catalyst Speciation Profile' } },
            scales: {
              x: { type: 'linear', min: 0, max: 100, title: { display: true, text: 'Catalyst ee (%)' } },
              y: { title: { display: true, text: 'Equilibrium Concentration (M)' } }
            }
          }
        });
      } else {
        chartSpecInstance.data.datasets = datasetsSpec;
        chartSpecInstance.update('none');
      }
    }

    function applyPreset(dg, khomo, ee0, g) {
      document.getElementById('slider_dg').value = dg;
      document.getElementById('slider_khomo').value = Math.log10(khomo);
      document.getElementById('slider_ee0').value = ee0;
      document.getElementById('slider_g').value = g;
      computeModel();
    }

    function handleCSV(event) {
      const file = event.target.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = function(e) {
        const lines = e.target.result.split('\\n');
        const pts = [];
        for (let line of lines) {
          line = line.trim();
          if (!line || line.startsWith('#')) continue;
          let parts = line.split(',');
          if (parts.length >= 2) {
            let x = parseFloat(parts[0]);
            let y = parseFloat(parts[1]);
            if (!isNaN(x) && !isNaN(y)) {
              if (x <= 1.0 && y <= 1.0 && (x > 0 || y > 0)) { x *= 100; y *= 100; }
              pts.push({ x: x, y: y });
            }
          }
        }
        experimentalData = pts;
        computeModel();
      };
      reader.readAsText(file);
    }

    ['slider_dg', 'slider_khomo', 'slider_ee0', 'slider_g', 'slider_temp'].forEach(id => {
      document.getElementById(id).addEventListener('input', computeModel);
    });

    window.onload = computeModel;
  </script>
</body>
</html>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    return output_path
