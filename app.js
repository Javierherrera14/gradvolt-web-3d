/**
 * ============================================================================
 * GradVolt Web 3D - Computational Physics Simulation Platform
 * Based on University Physics with Modern Physics by Sears & Zemansky
 * Chapters 21 (Coulomb & Electric Field), 22 (Gauss's Law), 23 (Potential)
 * ============================================================================
 */

(function () {
  'use strict';

  // --- Physical Constants (SI Units) ---
  const K_COULOMB = 8.9875517923e9; // N·m²/C² (k = 1 / (4πε₀))
  const EPSILON_0 = 8.8541878176e-12; // C²/(N·m²)
  const EPS_SMOOTH = 0.05; // Softening parameter (meters) to avoid singularities

  // --- Simulation State ---
  const state = {
    charges: [
      { id: 1, q: 5.0, x: -0.40, y: 0.00, z: 0.00 },
      { id: 2, q: -5.0, x: 0.40, y: 0.00, z: 0.00 }
    ],
    nextChargeId: 3,
    placingCharge: false,
    radius: 0.50, // Gaussian surface radius (m)
    resolution: 41, // N x N mesh grid
    activeTab: 'module-1',
    m1ViewMode: '3d', // '3d' or 'contour'
    m2ViewMode: '3d', // '3d' or '2d'
    m3ViewMode: '3d', // '3d' or '2d'
    selectedChargeA: 1,
    selectedChargeB: 2,
    colorscale: 'Electric',
    showEquipotential: true,
    showVectors: true,
    vectorDensity: 19,
    importedFromCsv: false,
    csvMesh: null
  };

  // --- Cached DOM Elements ---
  const elements = {
    // Charge management
    chargesCountBadge: document.getElementById('charges-count-badge'),
    btnCountDisplay: document.getElementById('btn-count-display'),
    btnAddCharge: document.getElementById('btn-add-charge'),
    chargesListContainer: document.getElementById('charges-list-container'),
    chargesPositionSummary: document.getElementById('charges-position-summary'),

    // Gauss Radius & resolution
    sliderRadius: document.getElementById('slider-radius'),
    radiusVal: document.getElementById('radius-val'),
    gaussRBadge: document.getElementById('gauss-r-badge'),
    selectResolution: document.getElementById('select-resolution'),
    labelResolutionStatus: document.getElementById('label-resolution-status'),
    btnResetDefaults: document.getElementById('btn-reset-defaults'),

    // Module 1
    btnView3D: document.getElementById('btn-view-3d'),
    btnViewContour: document.getElementById('btn-view-contour'),
    selectColorscale: document.getElementById('select-colorscale'),
    m1Vmax: document.getElementById('m1-vmax'),
    m1Vmin: document.getElementById('m1-vmin'),
    m1Vorigin: document.getElementById('m1-vorigin'),
    plotM1: document.getElementById('plot-module-1'),
    toggleTableBtn: document.getElementById('toggle-table-btn'),
    tableContainer: document.getElementById('table-container'),
    tableChevron: document.getElementById('table-chevron'),
    tablePointsCount: document.getElementById('table-points-count'),
    csvPreviewBody: document.getElementById('csv-preview-body'),

    // Module 2
    selectCoulombA: document.getElementById('select-coulomb-a'),
    selectCoulombB: document.getElementById('select-coulomb-b'),
    coulombDist: document.getElementById('coulomb-dist'),
    coulombForce: document.getElementById('coulomb-force'),
    coulombNature: document.getElementById('coulomb-nature'),
    coulombFnet: document.getElementById('coulomb-fnet'),
    coulombFnetComp: document.getElementById('coulomb-fnet-comp'),
    btnM2View3D: document.getElementById('btn-m2-view-3d'),
    btnM2View2D: document.getElementById('btn-m2-view-2d'),
    checkShowEquipotential: document.getElementById('check-show-equipotential'),
    checkShowVectors: document.getElementById('check-show-vectors'),
    selectVectorDensity: document.getElementById('select-vector-density'),
    plotM2: document.getElementById('plot-module-2'),

    // Module 3
    gaussInvarianceMsg: document.getElementById('gauss-invariance-msg'),
    gaussQenc: document.getElementById('gauss-qenc'),
    gaussQencSi: document.getElementById('gauss-qenc-si'),
    gaussFlux: document.getElementById('gauss-flux'),
    gaussArea: document.getElementById('gauss-area'),
    gaussRadiusDisp: document.getElementById('gauss-radius-disp'),
    gaussEavg: document.getElementById('gauss-eavg'),
    btnM3View3D: document.getElementById('btn-m3-view-3d'),
    btnM3View2D: document.getElementById('btn-m3-view-2d'),
    gaussChargesStatusList: document.getElementById('gauss-charges-status-list'),
    plotM3: document.getElementById('plot-module-3'),

    // CSV & Scenarios
    btnExportCsv: document.getElementById('btn-export-csv'),
    btnExportCsvSub: document.getElementById('btn-export-csv-sub'),
    inputImportCsv: document.getElementById('input-import-csv'),
    btnScenario1: document.getElementById('btn-scenario-1'),
    btnScenario2: document.getElementById('btn-scenario-2'),
    btnScenario3: document.getElementById('btn-scenario-3'),

    // Theory Modal
    btnOpenTheory: document.getElementById('btn-open-theory'),
    btnCloseTheory: document.getElementById('btn-close-theory'),
    btnCloseTheory2: document.getElementById('btn-close-theory-2'),
    modalTheory: document.getElementById('modal-theory'),
    toastContainer: document.getElementById('toast-container')
  };

  // --- Computational Physics Functions ---

  /**
   * Generates a uniform 1D array from min to max with N points
   */
  function linspace(min, max, N) {
    const arr = new Float64Array(N);
    const step = (max - min) / (N - 1);
    for (let i = 0; i < N; i++) {
      arr[i] = min + i * step;
    }
    return arr;
  }

  /**
   * Evaluates electric potential V at point (x, y, z) with superposition of N charges
   * V(x, y, z) = sum_i [ k * q_i / sqrt((x - xi)^2 + (y - yi)^2 + z^2 + eps^2) ]
   */
  function evaluatePotentialAt(x, y, z = 0) {
    let vTotal = 0;
    const len = state.charges.length;
    for (let i = 0; i < len; i++) {
      const c = state.charges[i];
      const qSI = c.q * 1e-9;
      const r = Math.sqrt((x - c.x) ** 2 + (y - c.y) ** 2 + z * z + EPS_SMOOTH ** 2);
      vTotal += K_COULOMB * (qSI / r);
    }
    return vTotal;
  }

  /**
   * Evaluates 3D Electric Field vector E at (x, y, z) analytically
   * E = sum_i [ k * q_i * (r - r_i) / |r - r_i|^3 ]
   */
  function evaluateElectricFieldAt(x, y, z = 0) {
    let ex = 0, ey = 0, ez = 0;
    const len = state.charges.length;
    for (let i = 0; i < len; i++) {
      const c = state.charges[i];
      const qSI = c.q * 1e-9;
      const dx = x - c.x;
      const dy = y - c.y;
      const dz = z - 0; // charges are always at z = 0
      const r2 = dx * dx + dy * dy + dz * dz + EPS_SMOOTH * EPS_SMOOTH;
      const r = Math.sqrt(r2);
      const factor = (K_COULOMB * qSI) / (r2 * r);
      ex += factor * dx;
      ey += factor * dy;
      ez += factor * dz;
    }
    const mag = Math.sqrt(ex * ex + ey * ey + ez * ez);
    return { ex, ey, ez, mag };
  }

  /**
   * Generates scalar potential matrix V[j][i] over 2D mesh at slice z=0
   */
  function computePotentialGrid(xArr, yArr) {
    const Ny = yArr.length;
    const Nx = xArr.length;
    const V = new Array(Ny);

    let vMin = Infinity;
    let vMax = -Infinity;

    for (let j = 0; j < Ny; j++) {
      const row = new Float64Array(Nx);
      const y = yArr[j];
      for (let i = 0; i < Nx; i++) {
        const x = xArr[i];
        const val = evaluatePotentialAt(x, y, 0);
        row[i] = val;
        if (val < vMin) vMin = val;
        if (val > vMax) vMax = val;
      }
      V[j] = row;
    }

    return { V, vMin, vMax };
  }

  /**
   * Formats numbers into scientific or readable human notation
   */
  function formatSci(val, unit = '', decimals = 2) {
    if (Math.abs(val) === 0 || isNaN(val)) return `0.00 ${unit}`.trim();
    const absVal = Math.abs(val);
    if (absVal >= 0.01 && absVal < 10000) {
      return `${val.toFixed(decimals)} ${unit}`.trim();
    }
    const exp = Math.floor(Math.log10(absVal));
    const mantissa = val / Math.pow(10, exp);
    return `${mantissa.toFixed(decimals)} × 10${toSuperscript(exp)} ${unit}`.trim();
  }

  function toSuperscript(num) {
    const digits = {
      '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³',
      '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'
    };
    return String(num).split('').map(c => digits[c] || c).join('');
  }

  // --- Plotly Theme & Layout Defaults ---
  const plotlyCommonLayout = {
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(15, 23, 42, 0.4)',
    font: {
      family: 'Inter, sans-serif',
      color: '#94a3b8'
    },
    margin: { t: 30, r: 30, b: 40, l: 50 },
    autosize: true
  };

  const plotlyConfig = {
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    displaylogo: false
  };

  // --- Rendering Pipeline ---
  let renderPending = false;

  function scheduleUpdate() {
    if (!renderPending) {
      renderPending = true;
      requestAnimationFrame(() => {
        updateSimulation();
        renderPending = false;
      });
    }
  }

  function updateSimulation() {
    updateUIElements();

    // 1. Grid setup
    const N = state.resolution;
    const xArr = linspace(-1.0, 1.0, N);
    const yArr = linspace(-1.0, 1.0, N);

    let V, vMin, vMax;
    if (state.importedFromCsv && state.csvMesh) {
      V = state.csvMesh.V;
      vMin = state.csvMesh.vMin;
      vMax = state.csvMesh.vMax;
    } else {
      const potResult = computePotentialGrid(xArr, yArr);
      V = potResult.V;
      vMin = potResult.vMin;
      vMax = potResult.vMax;
    }

    // 2. Metrics for Module 1
    const vOrigin = evaluatePotentialAt(0, 0, 0);
    if (elements.m1Vmax) elements.m1Vmax.textContent = `${vMax >= 0 ? '+' : ''}${formatSci(vMax, 'V')}`;
    if (elements.m1Vmin) elements.m1Vmin.textContent = `${formatSci(vMin, 'V')}`;
    if (elements.m1Vorigin) elements.m1Vorigin.textContent = `${formatSci(vOrigin, 'V')}`;

    // 3. Coulomb Calculation & Net Force (Module 2)
    updateCoulombMetrics();

    // 4. Gauss Calculation (Module 3)
    updateGaussMetrics();

    // 5. Update Active Module Visualization
    if (state.activeTab === 'module-1') {
      renderModule1Plot(xArr, yArr, V);
    } else if (state.activeTab === 'module-2') {
      renderModule2Plot(xArr, yArr, V);
    } else if (state.activeTab === 'module-3') {
      renderModule3Plot(xArr, yArr, V);
    }

    // 6. Update CSV Preview Table
    updatePreviewTable(xArr, yArr, V);
  }

  /**
   * Updates Coulomb Live Ribbon (Pairwise interaction + Net Force on Charge A)
   */
  function updateCoulombMetrics() {
    updateCoulombDropdowns();

    const cA = state.charges.find(c => c.id === state.selectedChargeA) || state.charges[0];
    const cB = state.charges.find(c => c.id === state.selectedChargeB) || (state.charges[1] || state.charges[0]);

    if (!cA) return;

    // Pairwise interaction between A and B
    if (cA && cB && cA.id !== cB.id) {
      const dx = cA.x - cB.x;
      const dy = cA.y - cB.y;
      const rAB = Math.sqrt(dx * dx + dy * dy);
      const qA_SI = cA.q * 1e-9;
      const qB_SI = cB.q * 1e-9;
      const forcePair = (K_COULOMB * Math.abs(qA_SI * qB_SI)) / Math.max(rAB * rAB, 1e-6);

      if (elements.coulombDist) elements.coulombDist.textContent = `${rAB.toFixed(3)} m`;
      if (elements.coulombForce) elements.coulombForce.textContent = formatSci(forcePair, 'N');

      if (elements.coulombNature) {
        if (cA.q === 0 || cB.q === 0) {
          elements.coulombNature.textContent = 'Nula';
          elements.coulombNature.className = 'text-[10px] font-semibold text-slate-400';
        } else if (cA.q * cB.q > 0) {
          elements.coulombNature.textContent = 'Repulsiva';
          elements.coulombNature.className = 'text-[10px] font-semibold text-rose-400';
        } else {
          elements.coulombNature.textContent = 'Atractiva';
          elements.coulombNature.className = 'text-[10px] font-semibold text-cyan-400';
        }
      }
    } else {
      if (elements.coulombDist) elements.coulombDist.textContent = '0.000 m';
      if (elements.coulombForce) elements.coulombForce.textContent = '0.00 N';
      if (elements.coulombNature) {
        elements.coulombNature.textContent = 'Misma carga';
        elements.coulombNature.className = 'text-[10px] font-semibold text-slate-500';
      }
    }

    // Net Force on Charge A due to ALL other charges
    let fNetX = 0, fNetY = 0;
    const qA_SI = cA.q * 1e-9;

    state.charges.forEach(cj => {
      if (cj.id !== cA.id) {
        const dx = cA.x - cj.x;
        const dy = cA.y - cj.y;
        const r2 = dx * dx + dy * dy;
        const r = Math.sqrt(r2);
        if (r > 1e-6) {
          const qj_SI = cj.q * 1e-9;
          const factor = (K_COULOMB * qA_SI * qj_SI) / (r2 * r);
          fNetX += factor * dx;
          fNetY += factor * dy;
        }
      }
    });

    const fNetMag = Math.sqrt(fNetX * fNetX + fNetY * fNetY);
    if (elements.coulombFnet) elements.coulombFnet.textContent = formatSci(fNetMag, 'N');
    if (elements.coulombFnetComp) {
      elements.coulombFnetComp.textContent = `(${formatSci(fNetX, '', 1)}, ${formatSci(fNetY, '', 1)}) N`;
    }
  }

  /**
   * Updates Gauss Metrics & Badges
   */
  function updateGaussMetrics() {
    let qEnc = 0;
    const R = state.radius;

    // Badges HTML
    let badgesHtml = '';

    state.charges.forEach((c, idx) => {
      const dist = Math.sqrt(c.x * c.x + c.y * c.y);
      const isInside = dist <= R;

      if (isInside) {
        qEnc += c.q;
      }

      const qColor = c.q > 0 ? 'text-rose-300' : (c.q < 0 ? 'text-cyan-300' : 'text-slate-400');
      const badgeBorder = isInside
        ? 'bg-emerald-500/20 text-emerald-300 border-emerald-500/40 shadow-[0_0_8px_rgba(16,185,129,0.25)]'
        : 'bg-slate-800 text-slate-500 border-slate-700';

      badgesHtml += `
        <div class="px-2.5 py-1 rounded-full text-[11px] font-semibold border flex items-center gap-1.5 ${badgeBorder}">
          <span class="w-2 h-2 rounded-full ${isInside ? 'bg-emerald-400' : 'bg-slate-600'}"></span>
          <span class="${qColor}">q${idx + 1} (${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC):</span>
          <span class="font-bold">${isInside ? 'DENTRO' : 'FUERA'}</span>
        </div>
      `;
    });

    if (elements.gaussChargesStatusList) {
      elements.gaussChargesStatusList.innerHTML = badgesHtml;
    }

    const qEncSI = qEnc * 1e-9;
    const fluxTheoretical = qEncSI / EPSILON_0;
    const areaSphere = 4 * Math.PI * (R * R);
    const eAvg = Math.abs(fluxTheoretical) / areaSphere;

    if (elements.gaussQenc) elements.gaussQenc.textContent = `${qEnc >= 0 ? '+' : ''}${qEnc.toFixed(2)} nC`;
    if (elements.gaussQencSi) elements.gaussQencSi.textContent = formatSci(qEncSI, 'C');
    if (elements.gaussFlux) elements.gaussFlux.textContent = formatSci(fluxTheoretical, 'N·m²/C');
    if (elements.gaussArea) elements.gaussArea.textContent = `${areaSphere.toFixed(3)} m²`;
    if (elements.gaussRadiusDisp) elements.gaussRadiusDisp.textContent = `R = ${R.toFixed(2)} m`;
    if (elements.gaussEavg) elements.gaussEavg.textContent = formatSci(eAvg, 'N/C');

    if (elements.gaussInvarianceMsg) {
      if (Math.abs(qEnc) < 1e-6) {
        elements.gaussInvarianceMsg.innerHTML = 'Carga neta interior <strong class="text-amber-300">Q_enc = 0</strong>: El flujo eléctrico neto es rigurosamente cero (<strong class="text-white">Φ_E = 0</strong>).';
      } else {
        elements.gaussInvarianceMsg.innerHTML = `Flujo constante en <strong class="text-emerald-300">${formatSci(fluxTheoretical, 'N·m²/C')}</strong>. Al mover el radio R, el flujo permanece idéntico mientras ninguna carga cruce la frontera.`;
      }
    }
  }

  /**
   * Syncs Global UI Elements with State
   */
  function updateUIElements() {
    if (elements.chargesCountBadge) {
      elements.chargesCountBadge.textContent = `${state.charges.length} / 10`;
    }
    if (elements.btnCountDisplay) {
      elements.btnCountDisplay.textContent = state.charges.length;
    }

    // Gauss Radius
    if (elements.sliderRadius) elements.sliderRadius.value = state.radius;
    if (elements.radiusVal) elements.radiusVal.textContent = `${state.radius.toFixed(2)} m`;
    if (elements.gaussRBadge) elements.gaussRBadge.textContent = `R = ${state.radius.toFixed(2)} m`;

    // Colorscale & Resolution
    if (elements.selectColorscale) elements.selectColorscale.value = state.colorscale;
    if (elements.labelResolutionStatus) elements.labelResolutionStatus.textContent = `Malla ${state.resolution}×${state.resolution}`;
     updatePositionSummary();
}

/**
 * Actualiza el panel superior con los vectores de posición de cada carga
 */
function updatePositionSummary() {
  const container = elements.chargesPositionSummary;
  if (!container) return;

  let html = '';
  state.charges.forEach((c, idx) => {
    const isPos = c.q > 0;
    const isNeg = c.q < 0;
    const borderClass = isPos
      ? 'border-rose-500/40 bg-rose-950/20 text-rose-200'
      : (isNeg ? 'border-cyan-500/40 bg-cyan-950/20 text-cyan-200' : 'border-slate-700 bg-slate-900 text-slate-300');
    const qColor = isPos ? 'text-rose-400' : (isNeg ? 'text-cyan-400' : 'text-slate-400');

    html += `
      <div class="px-2 py-1 rounded-lg border ${borderClass} text-[10px] font-mono flex items-center justify-between">
        <span class="font-bold ${qColor}">r_${idx + 1}:</span>
        <span>(${c.x >= 0 ? '+' : ''}${c.x.toFixed(2)}, ${c.y >= 0 ? '+' : ''}${c.y.toFixed(2)}) m</span>
      </div>
    `;
  });

  container.innerHTML = html;
}

  /**
   * Renders Dynamic Charge Cards in the Sidebar (Single List, q, x, y Sliders)
   */
  function renderChargesList() {
    const container = elements.chargesListContainer;
    if (!container) return;

    let html = '';

    state.charges.forEach((c, idx) => {
      const chargeNum = idx + 1;
      const isPositive = c.q > 0;
      const isNegative = c.q < 0;
      const cardType = isPositive ? 'charge-positive' : (isNegative ? 'charge-negative' : 'charge-neutral');
      const dotColor = isPositive ? 'bg-rose-500 shadow-[0_0_8px_rgba(244,63,94,0.8)]' : (isNegative ? 'bg-cyan-400 shadow-[0_0_8px_rgba(6,182,212,0.8)]' : 'bg-slate-400');
      const badgeBg = isPositive ? 'bg-rose-500/20 text-rose-300 border-rose-500/30' : (isNegative ? 'bg-cyan-500/20 text-cyan-300 border-cyan-500/30' : 'bg-slate-800 text-slate-400 border-slate-700');
      const sliderAccent = isPositive ? 'slider-accent-pos' : (isNegative ? 'slider-accent-neg' : 'slider-accent-neutral');

      html += `
        <div class="glass-card charge-card ${cardType} rounded-2xl p-3.5 border relative" data-charge-id="${c.id}">
          <div class="flex items-center justify-between mb-2.5">
            <div class="flex items-center gap-2">
              <span class="w-3 h-3 rounded-full ${dotColor}"></span>
              <h3 class="text-xs font-bold tracking-wide text-slate-200">
                Carga Puntual $q_{${chargeNum}}$
              </h3>
            </div>
            
            <div class="flex items-center gap-1.5">
              <span class="charge-badge-text text-[11px] font-mono font-bold px-2 py-0.5 rounded ${badgeBg} border">
                ${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC
              </span>
              ${state.charges.length > 1 ? `
                <button class="btn-delete-charge p-1 text-slate-500 hover:text-rose-400 transition-colors" data-id="${c.id}" title="Eliminar carga q${chargeNum}">
                  <i class="fa-solid fa-trash-can text-xs"></i>
                </button>
              ` : ''}
            </div>
          </div>

          <div class="space-y-2.5">
            <!-- Magnitude slider -->
            <div>
              <div class="flex justify-between text-[11px] text-slate-400 mb-0.5">
                <span>Magnitud ($q_{${chargeNum}}$)</span>
                <span class="font-mono text-slate-200 font-bold charge-label-q">${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC</span>
              </div>
              <input type="range" min="-10" max="10" step="0.5" value="${c.q}" class="w-full range-slider ${sliderAccent} charge-slider-q" data-id="${c.id}">
              <div class="flex justify-between text-[9px] text-slate-500 font-mono mt-0.5">
                <span>-10 nC</span>
                <span>0</span>
                <span>+10 nC</span>
              </div>
            </div>

            <!-- Position Sliders (x, y) -->
            <div class="grid grid-cols-2 gap-2 pt-1 border-t border-slate-800/80">
              <div>
                <div class="flex justify-between text-[10px] text-slate-400 mb-0.5">
                  <span>Pos $x_{${chargeNum}}$</span>
                  <span class="font-mono text-slate-300 charge-label-x">${c.x.toFixed(2)} m</span>
                </div>
                <input type="range" min="-0.9" max="0.9" step="0.05" value="${c.x}" class="w-full range-slider ${sliderAccent} charge-slider-x" data-id="${c.id}">
              </div>

              <div>
                <div class="flex justify-between text-[10px] text-slate-400 mb-0.5">
                  <span>Pos $y_{${chargeNum}}$</span>
                  <span class="font-mono text-slate-300 charge-label-y">${c.y.toFixed(2)} m</span>
                </div>
                <input type="range" min="-0.9" max="0.9" step="0.05" value="${c.y}" class="w-full range-slider ${sliderAccent} charge-slider-y" data-id="${c.id}">
              </div>
            </div>
          </div>
        </div>
      `;
    });

    container.innerHTML = html;

    // Attach listeners to dynamically generated sliders
    container.querySelectorAll('.charge-slider-q').forEach(input => {
      input.addEventListener('input', e => {
        const id = parseInt(e.target.getAttribute('data-id'), 10);
        const c = state.charges.find(ch => ch.id === id);
        if (c) {
          c.q = parseFloat(e.target.value);
          state.importedFromCsv = false;

          // Smooth live label update without redrawing whole card
          const card = e.target.closest('.charge-card');
          if (card) {
            const labelQ = card.querySelector('.charge-label-q');
            if (labelQ) labelQ.textContent = `${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC`;
            const badge = card.querySelector('.charge-badge-text');
            if (badge) badge.textContent = `${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC`;
          }
          scheduleUpdate();
        }
      });

      input.addEventListener('change', () => {
        // Full refresh on release to update colors and styling if sign changed
        renderChargesList();
      });
    });

    container.querySelectorAll('.charge-slider-x').forEach(input => {
      input.addEventListener('input', e => {
        const id = parseInt(e.target.getAttribute('data-id'), 10);
        const c = state.charges.find(ch => ch.id === id);
        if (c) {
          c.x = parseFloat(e.target.value);
          state.importedFromCsv = false;
          const card = e.target.closest('.charge-card');
          if (card) {
            const labelX = card.querySelector('.charge-label-x');
            if (labelX) labelX.textContent = `${c.x.toFixed(2)} m`;
          }
          scheduleUpdate();
        }
      });
    });

    container.querySelectorAll('.charge-slider-y').forEach(input => {
      input.addEventListener('input', e => {
        const id = parseInt(e.target.getAttribute('data-id'), 10);
        const c = state.charges.find(ch => ch.id === id);
        if (c) {
          c.y = parseFloat(e.target.value);
          state.importedFromCsv = false;
          const card = e.target.closest('.charge-card');
          if (card) {
            const labelY = card.querySelector('.charge-label-y');
            if (labelY) labelY.textContent = `${c.y.toFixed(2)} m`;
          }
          scheduleUpdate();
        }
      });
    });

    // Delete buttons
    container.querySelectorAll('.btn-delete-charge').forEach(btn => {
      btn.addEventListener('click', () => {
        const id = parseInt(btn.getAttribute('data-id'), 10);
        deleteCharge(id);
      });
    });

    renderMathInElement();
  }

  /**
   * Updates Coulomb Pair Selectors
   */
  function updateCoulombDropdowns() {
    const selA = elements.selectCoulombA;
    const selB = elements.selectCoulombB;
    if (!selA || !selB) return;

    const currentA = state.selectedChargeA;
    const currentB = state.selectedChargeB;

    let optionsHtml = '';
    state.charges.forEach((c, idx) => {
      optionsHtml += `<option value="${c.id}">q${idx + 1} (${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC)</option>`;
    });

    selA.innerHTML = optionsHtml;
    selB.innerHTML = optionsHtml;

    if (state.charges.some(c => c.id === currentA)) {
      selA.value = String(currentA);
    } else {
      selA.value = String(state.charges[0].id);
      state.selectedChargeA = state.charges[0].id;
    }

    if (state.charges.some(c => c.id === currentB)) {
      selB.value = String(currentB);
    } else {
      const defaultB = state.charges.length > 1 ? state.charges[1].id : state.charges[0].id;
      selB.value = String(defaultB);
      state.selectedChargeB = defaultB;
    }
  }

  /**
   * Updates Add Charge button UI and cursor according to placingCharge state
   */
  function updateAddButtonUI() {
    const btn = elements.btnAddCharge;
    const plots = [elements.plotM1, elements.plotM2, elements.plotM3];

    if (state.placingCharge) {
      if (btn) {
        btn.innerHTML = `<span>📍 Haz clic en el mapa para ubicar (o clic aquí para pos. automática)</span>`;
        btn.classList.add('ring-2', 'ring-brand-400');
      }
      plots.forEach(el => {
        if (el) {
          el.style.cursor = 'crosshair';
          el.classList.add('placing-charge-active');
        }
      });
    } else {
      if (btn) {
        btn.innerHTML = `<i class="fa-solid fa-circle-plus mr-1"></i><span>+ Agregar Carga (<span id="btn-count-display">${state.charges.length}</span>/10)</span>`;
        btn.classList.remove('ring-2', 'ring-brand-400');
        elements.btnCountDisplay = document.getElementById('btn-count-display');
      }
      plots.forEach(el => {
        if (el) {
          el.style.cursor = '';
          el.classList.remove('placing-charge-active');
        }
      });
    }
  }

  /**
   * Attaches Plotly click handler to plot element to place charges interactively
   * Note: Invoked strictly after Plotly.newPlot to ensure plotEl is initialized.
   */
  function attachPlotClickHandler(plotEl) {
    if (!plotEl || typeof plotEl.on !== 'function') return;

    if (typeof plotEl.removeAllListeners === 'function') {
      plotEl.removeAllListeners("plotly_click");
    }

    plotEl.on("plotly_click", function (data) {
      if (!state.placingCharge) return;
      if (!data || !data.points || !data.points.length) return;

      const pt = data.points[0];
      if (typeof pt.x !== 'number' || typeof pt.y !== 'number') return;
      if (isNaN(pt.x) || isNaN(pt.y)) return;

      const clampAndRound = val => Math.max(-0.90, Math.min(0.90, Math.round(val / 0.05) * 0.05));
      const newX = parseFloat(clampAndRound(pt.x).toFixed(2));
      const newY = parseFloat(clampAndRound(pt.y).toFixed(2));
      const newQ = (state.charges.length % 2 === 0) ? 5.0 : -5.0;

      const added = addNewCharge(newQ, newX, newY);
      state.placingCharge = false;
      updateAddButtonUI();

      if (added) {
        showToast(`Carga puntual q${state.charges.length} ubicada en (${newX >= 0 ? '+' : ''}${newX}, ${newY >= 0 ? '+' : ''}${newY}) m.`, 'success');
      }
    });
  }

  /**
   * Adds a new charge with specific parameters
   */
  function addNewCharge(q, x, y) {
    if (state.charges.length >= 10) {
      showToast('Límite máximo de 10 cargas alcanzado.', 'error');
      return false;
    }

    const newCharge = {
      id: state.nextChargeId++,
      q: parseFloat(q),
      x: parseFloat(x),
      y: parseFloat(y),
      z: 0.00
    };

    state.charges.push(newCharge);
    state.importedFromCsv = false;
    renderChargesList();
    updateCoulombDropdowns();
    updateAddButtonUI();
    scheduleUpdate();
    return true;
  }

  /**
   * Direct handler for the Add Charge button
   */
  function handleAddChargeClick() {
    if (state.charges.length >= 10) {
      showToast('Límite máximo de 10 cargas alcanzado.', 'error');
      return;
    }

    if (!state.placingCharge) {
      state.placingCharge = true;
      updateAddButtonUI();
      showToast('Haz clic sobre la gráfica para colocar la carga puntual (o clic aquí para posición automática).', 'info');
    } else {
      state.placingCharge = false;
      const count = state.charges.length;
      // Distribute nicely along an orbit to avoid placing right on top of another charge
      const angle = (count * 2 * Math.PI) / 8;
      const newX = parseFloat((0.45 * Math.cos(angle)).toFixed(2));
      const newY = parseFloat((0.45 * Math.sin(angle)).toFixed(2));
      const newQ = (count % 2 === 0) ? 5.0 : -5.0;

      const added = addNewCharge(newQ, newX, newY);
      updateAddButtonUI();
      if (added) {
        showToast(`Carga puntual q${state.charges.length} agregada (${newQ >= 0 ? '+' : ''}${newQ} nC en x=${newX}, y=${newY} m).`, 'success');
      }
    }
  }

  /**
   * Deletes a charge by ID
   */
  function deleteCharge(id) {
    if (state.charges.length <= 1) {
      showToast('Debe permanecer al menos 1 carga activa.', 'error');
      return;
    }

    state.charges = state.charges.filter(c => c.id !== id);
    if (state.selectedChargeA === id) state.selectedChargeA = state.charges[0].id;
    if (state.selectedChargeB === id) state.selectedChargeB = state.charges.length > 1 ? state.charges[1].id : state.charges[0].id;

    renderChargesList();
    updateCoulombDropdowns();
    updateAddButtonUI();
    scheduleUpdate();
    showToast('Carga eliminada correctamente.', 'info');
  }

  // --- Module 1 Plot: 3D Surface / 2D Contours (Robust 3D <-> 2D Purge Transition) ---
  function renderModule1Plot(xArr, yArr, V) {
    const plotEl = elements.plotM1;
    if (!plotEl) return;

    const targetMode = state.m1ViewMode; // '3d' or 'contour'
    const data = [];

    if (targetMode === '3d') {
      // 3D Surface
      data.push({
        type: 'surface',
        x: Array.from(xArr),
        y: Array.from(yArr),
        z: V.map(row => Array.from(row)),
        colorscale: state.colorscale,
        contours: {
          z: {
            show: true,
            usecolormap: true,
            highlightcolor: '#ffffff',
            project: { z: true }
          }
        },
        lighting: {
          roughness: 0.3,
          ambient: 0.8,
          diffuse: 0.9,
          specular: 0.4
        },
        colorbar: {
          title: { text: 'V (Voltios)', side: 'right', font: { color: '#e2e8f0', size: 12 } },
          tickfont: { color: '#94a3b8', size: 10 },
          len: 0.75,
          thickness: 16
        }
      });

      // Scatter 3D markers for all charges projected onto surface
      data.push({
        type: 'scatter3d',
        mode: 'markers+text',
        x: state.charges.map(c => c.x),
        y: state.charges.map(c => c.y),
        z: state.charges.map(c => evaluatePotentialAt(c.x, c.y, 0)),
        text: state.charges.map((c, idx) => `q${idx + 1} (${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC)`),
        textposition: 'top center',
        textfont: { family: 'Inter', size: 11, color: '#ffffff' },
        marker: {
          size: 7,
          color: state.charges.map(c => c.q > 0 ? '#f43f5e' : (c.q < 0 ? '#06b6d4' : '#94a3b8')),
          line: { color: '#ffffff', width: 2 }
        },
        showlegend: false
      });

      const layout = {
        ...plotlyCommonLayout,
        scene: {
          dragmode: 'orbit',
          xaxis: { title: 'x (m)', color: '#94a3b8', gridcolor: '#334155', backgroundcolor: 'rgba(15, 23, 42, 0.6)' },
          yaxis: { title: 'y (m)', color: '#94a3b8', gridcolor: '#334155', backgroundcolor: 'rgba(15, 23, 42, 0.6)' },
          zaxis: { title: 'Potencial V (V)', color: '#94a3b8', gridcolor: '#334155', backgroundcolor: 'rgba(15, 23, 42, 0.6)' },
          camera: { eye: { x: 1.45, y: -1.45, z: 1.25 } }
        }
      };

      if (plotEl._currentMode !== '3d') {
        Plotly.purge(plotEl);
        Plotly.newPlot(plotEl, data, layout, plotlyConfig);
        plotEl._currentMode = '3d';
        attachPlotClickHandler(plotEl);
      } else {
        Plotly.react(plotEl, data, layout, plotlyConfig);
      }
    } else {
      // 2D Contour (Equipotential lines)
      data.push({
        type: 'contour',
        x: Array.from(xArr),
        y: Array.from(yArr),
        z: V.map(row => Array.from(row)),
        colorscale: state.colorscale,
        ncontours: 28,
        contours: {
          coloring: 'heatmap',
          showlabels: true,
          labelfont: { family: 'Inter', size: 10, color: '#ffffff' }
        },
        line: { width: 0.8, color: 'rgba(255,255,255,0.4)' },
        colorbar: {
          title: { text: 'V (Voltios)', font: { color: '#e2e8f0', size: 12 } },
          tickfont: { color: '#94a3b8', size: 10 },
          thickness: 16
        }
      });

      // Overlay all active charges in 2D
      data.push({
        type: 'scatter',
        mode: 'markers+text',
        x: state.charges.map(c => c.x),
        y: state.charges.map(c => c.y),
        text: state.charges.map((c, idx) => `q${idx + 1} (${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC)`),
        textposition: 'top center',
        textfont: { family: 'Inter', size: 11, color: '#f8fafc' },
        marker: {
          size: 16,
          color: state.charges.map(c => c.q > 0 ? '#f43f5e' : (c.q < 0 ? '#06b6d4' : '#94a3b8')),
          line: { color: '#ffffff', width: 2 }
        },
        showlegend: false
      });

      const layout = {
        ...plotlyCommonLayout,
        xaxis: { title: 'x (metros)', range: [-1.05, 1.05], zerolinecolor: '#475569', gridcolor: '#1e293b' },
        yaxis: { title: 'y (metros)', range: [-1.05, 1.05], zerolinecolor: '#475569', gridcolor: '#1e293b', scaleanchor: 'x', scaleratio: 1 }
      };

      if (plotEl._currentMode !== 'contour') {
        Plotly.purge(plotEl);
        Plotly.newPlot(plotEl, data, layout, plotlyConfig);
        plotEl._currentMode = 'contour';
        attachPlotClickHandler(plotEl);
      } else {
        Plotly.react(plotEl, data, layout, plotlyConfig);
      }
    }
  }

  /**
   * Generates continuous electric field streamlines via step-by-step RK2 integration (ds = 0.025)
   * with directional arrowheads at 35% and 70% of each curve indicating the exact E field direction.
   * Sears & Zemansky Chapters 21-23.
   */
  function buildFieldStreamlines() {
    const ds = 0.025; // integration step in meters
    const maxSteps = 120; // maximum RK2 steps per streamline
    const linesPerCharge = 12;
    const launchRadius = 0.045; // initial offset from charge center in meters
    const collisionDist = 0.06; // stop when within 0.06 m of another charge

    const posCharges = state.charges.filter(c => c.q > 0.01);
    const negCharges = state.charges.filter(c => c.q < -0.01);

    // If positive charges exist, emanate forward from positive charges.
    // If only negative charges exist, emanate backwards from negative charges.
    let seedCharges = [];
    let isReverse = false;

    if (posCharges.length > 0) {
      seedCharges = posCharges;
      isReverse = false;
    } else if (negCharges.length > 0) {
      seedCharges = negCharges;
      isReverse = true;
    } else {
      return { lineX: [], lineY: [], lineZ: [] };
    }

    const lineX = [];
    const lineY = [];
    const lineZ = [];

    // Evaluates normalized field direction
    function getUnitE(x, y, sign) {
      const field = evaluateElectricFieldAt(x, y, 0);
      if (field.mag < 1e-6) return null;
      return {
        ux: (field.ex / field.mag) * sign,
        uy: (field.ey / field.mag) * sign
      };
    }

    // Appends a directional arrowhead pointing along (ux, uy) with tip at (tipX, tipY)
    function addArrowHead(tipX, tipY, ux, uy) {
      const arrowLen = 0.032;
      const angle = 0.45;
      const w1x = tipX - arrowLen * (ux * Math.cos(angle) - uy * Math.sin(angle));
      const w1y = tipY - arrowLen * (uy * Math.cos(angle) + ux * Math.sin(angle));
      const w2x = tipX - arrowLen * (ux * Math.cos(-angle) - uy * Math.sin(-angle));
      const w2y = tipY - arrowLen * (uy * Math.cos(-angle) + ux * Math.sin(-angle));

      lineX.push(w1x, tipX, w2x, null);
      lineY.push(w1y, tipY, w2y, null);
      lineZ.push(0.0, 0.0, 0.0, null);
    }

    const sign = isReverse ? -1 : 1;

    seedCharges.forEach(seedC => {
      for (let k = 0; k < linesPerCharge; k++) {
        const theta = (k / linesPerCharge) * 2 * Math.PI;
        const x0 = seedC.x + launchRadius * Math.cos(theta);
        const y0 = seedC.y + launchRadius * Math.sin(theta);

        const curvePts = [{ x: x0, y: y0 }];
        let currX = x0;
        let currY = y0;

        for (let step = 0; step < maxSteps; step++) {
          // RK2 - Step 1: evaluate at current point
          const u1 = getUnitE(currX, currY, sign);
          if (!u1) break;

          // Midpoint evaluation
          const midX = currX + 0.5 * ds * u1.ux;
          const midY = currY + 0.5 * ds * u1.uy;

          // RK2 - Step 2: evaluate at midpoint
          const u2 = getUnitE(midX, midY, sign);
          if (!u2) break;

          const nextX = currX + ds * u2.ux;
          const nextY = currY + ds * u2.uy;

          // Check domain bounds [-1.0, 1.0]
          if (nextX < -1.0 || nextX > 1.0 || nextY < -1.0 || nextY > 1.0) {
            curvePts.push({ x: nextX, y: nextY });
            break;
          }

          // Check collision with another charge (< 0.06 m)
          let collided = false;
          for (const ch of state.charges) {
            if (ch.id === seedC.id && step < 2) continue; // ignore origin charge on takeoff
            const dist = Math.hypot(nextX - ch.x, nextY - ch.y);
            if (dist < collisionDist) {
              collided = true;
              break;
            }
          }

          curvePts.push({ x: nextX, y: nextY });
          if (collided) break;

          currX = nextX;
          currY = nextY;
        }

        if (curvePts.length < 2) continue;

        // If reverse integration from negative charges, reverse order so points flow toward charge
        if (isReverse) {
          curvePts.reverse();
        }

        // Add streamline points to line array
        for (let i = 0; i < curvePts.length; i++) {
          lineX.push(curvePts[i].x);
          lineY.push(curvePts[i].y);
          lineZ.push(0.0);
        }
        lineX.push(null);
        lineY.push(null);
        lineZ.push(null);

        // Add directional arrowheads at 35% and 70% of the streamline curve
        const nPts = curvePts.length;
        if (nPts >= 4) {
          const arrowIndices = [
            Math.max(1, Math.floor(nPts * 0.35)),
            Math.min(nPts - 2, Math.floor(nPts * 0.70))
          ];

          const uniqueIndices = Array.from(new Set(arrowIndices));

          uniqueIndices.forEach(idx => {
            const pt = curvePts[idx];
            const field = evaluateElectricFieldAt(pt.x, pt.y, 0);
            if (field.mag > 1e-6) {
              const ux = field.ex / field.mag;
              const uy = field.ey / field.mag;
              addArrowHead(pt.x, pt.y, ux, uy);
            }
          });
        }
      }
    });

    return { lineX, lineY, lineZ };
  }

  // --- Helper to build 2D vector quiver lines ---
  function buildQuiverVectors(density, xRange, yRange) {
    const qX = linspace(xRange[0], xRange[1], density);
    const qY = linspace(yRange[0], yRange[1], density);
    const arrowLen = 1.3 / density;
    const headLen = arrowLen * 0.35;
    const headAngle = 0.45;

    const lineX = [];
    const lineY = [];

    for (let j = 0; j < density; j++) {
      const y0 = qY[j];
      for (let i = 0; i < density; i++) {
        const x0 = qX[i];
        const { ex, ey, mag } = evaluateElectricFieldAt(x0, y0, 0);

        if (mag < 1e-4) continue;

        const ux = ex / mag;
        const uy = ey / mag;

        const startX = x0 - (ux * arrowLen * 0.4);
        const startY = y0 - (uy * arrowLen * 0.4);
        const endX = x0 + (ux * arrowLen * 0.6);
        const endY = y0 + (uy * arrowLen * 0.6);

        const b1X = endX - headLen * (ux * Math.cos(headAngle) - uy * Math.sin(headAngle));
        const b1Y = endY - headLen * (uy * Math.cos(headAngle) + ux * Math.sin(headAngle));
        const b2X = endX - headLen * (ux * Math.cos(-headAngle) - uy * Math.sin(-headAngle));
        const b2Y = endY - headLen * (uy * Math.cos(-headAngle) + ux * Math.sin(-headAngle));

        lineX.push(startX, endX, b1X, endX, b2X, null);
        lineY.push(startY, endY, b1Y, endY, b2Y, null);
      }
    }

    return { lineX, lineY };
  }

  // --- Module 2 Plot: 3D Rotatable Cone Field / 2D Quiver Field (Robust 3D <-> 2D Purge Transition) ---
  function renderModule2Plot(xArr, yArr, V) {
    const plotEl = elements.plotM2;
    if (!plotEl) return;

    const targetMode = state.m2ViewMode; // '3d' or '2d'
    const data = [];

    if (targetMode === '3d') {
      // 3D Rotatable Vector Field via Cones
      const coneDensity = 7; // 7x7x7 = 343 sample points in 3D
      const cRange = linspace(-0.80, 0.80, coneDensity);

      const cx = [], cy = [], cz = [];
      const cu = [], cv = [], cw = [];

      for (let iz = 0; iz < coneDensity; iz++) {
        const z = cRange[iz];
        for (let iy = 0; iy < coneDensity; iy++) {
          const y = cRange[iy];
          for (let ix = 0; ix < coneDensity; ix++) {
            const x = cRange[ix];

            const { ex, ey, ez, mag } = evaluateElectricFieldAt(x, y, z);
            if (mag > 1e-3) {
              cx.push(x);
              cy.push(y);
              cz.push(z);
              const scale = Math.min(1.0, Math.log10(mag + 1) * 0.25 + 0.3);
              cu.push((ex / mag) * scale);
              cv.push((ey / mag) * scale);
              cw.push((ez / mag) * scale);
            }
          }
        }
      }

      data.push({
        type: 'cone',
        x: cx,
        y: cy,
        z: cz,
        u: cu,
        v: cv,
        w: cw,
        sizemode: 'scaled',
        sizeref: 0.35,
        anchor: 'tail',
        colorscale: 'Electric',
        showscale: false
      });

      // 3D Continuous Field Streamlines with directional arrows at z = 0
      const streamlines3D = buildFieldStreamlines();
      if (streamlines3D.lineX.length > 0) {
        data.push({
          type: 'scatter3d',
          mode: 'lines',
          x: streamlines3D.lineX,
          y: streamlines3D.lineY,
          z: streamlines3D.lineZ,
          line: { color: '#38bdf8', width: 3 },
          name: 'Líneas de Campo Continuas',
          hoverinfo: 'none',
          showlegend: true
        });
      }

      // 3D Charge spheres at z = 0
      data.push({
        type: 'scatter3d',
        mode: 'markers+text',
        x: state.charges.map(c => c.x),
        y: state.charges.map(c => c.y),
        z: state.charges.map(() => 0.0),
        text: state.charges.map((c, idx) => `q${idx + 1} (${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC)`),
        textposition: 'top center',
        textfont: { family: 'Inter', size: 11, color: '#ffffff' },
        marker: {
          size: 8,
          color: state.charges.map(c => c.q > 0 ? '#f43f5e' : (c.q < 0 ? '#06b6d4' : '#94a3b8')),
          line: { color: '#ffffff', width: 2 }
        },
        showlegend: false
      });

      // Line between selected Charge A and B
      const cA = state.charges.find(c => c.id === state.selectedChargeA) || state.charges[0];
      const cB = state.charges.find(c => c.id === state.selectedChargeB) || (state.charges[1] || state.charges[0]);
      if (cA && cB && cA.id !== cB.id) {
        data.push({
          type: 'scatter3d',
          mode: 'lines',
          x: [cA.x, cB.x],
          y: [cA.y, cB.y],
          z: [0, 0],
          line: { color: '#f59e0b', width: 4, dash: 'dot' },
          name: 'Distancia r_AB',
          hoverinfo: 'none'
        });
      }

      const layout = {
        ...plotlyCommonLayout,
        scene: {
          dragmode: 'orbit',
          xaxis: { title: 'x (m)', color: '#94a3b8', gridcolor: '#334155', backgroundcolor: 'rgba(15, 23, 42, 0.6)' },
          yaxis: { title: 'y (m)', color: '#94a3b8', gridcolor: '#334155', backgroundcolor: 'rgba(15, 23, 42, 0.6)' },
          zaxis: { title: 'z (m)', color: '#94a3b8', gridcolor: '#334155', backgroundcolor: 'rgba(15, 23, 42, 0.6)' },
          camera: { eye: { x: 1.5, y: -1.5, z: 1.2 } }
        }
      };

      if (plotEl._currentMode !== '3d') {
        Plotly.purge(plotEl);
        Plotly.newPlot(plotEl, data, layout, plotlyConfig);
        plotEl._currentMode = '3d';
        attachPlotClickHandler(plotEl);
      } else {
        Plotly.react(plotEl, data, layout, plotlyConfig);
      }
    } else {
      // 2D Quiver Field View
      if (state.showEquipotential) {
        data.push({
          type: 'contour',
          x: Array.from(xArr),
          y: Array.from(yArr),
          z: V.map(row => Array.from(row)),
          colorscale: 'Viridis',
          opacity: 0.6,
          ncontours: 26,
          contours: {
            coloring: 'heatmap',
            showlabels: true,
            labelfont: { family: 'Inter', size: 9, color: 'rgba(255,255,255,0.85)' }
          },
          line: { width: 1, color: 'rgba(255,255,255,0.45)' },
          colorbar: { title: { text: 'V (V)', font: { color: '#e2e8f0', size: 11 } }, thickness: 14, len: 0.8 },
          hoverinfo: 'x+y+z'
        });
      }

      // Continuous Field Streamlines with directional arrows in 2D
      const streamlines2D = buildFieldStreamlines();
      if (streamlines2D.lineX.length > 0) {
        data.push({
          type: 'scatter',
          mode: 'lines',
          x: streamlines2D.lineX,
          y: streamlines2D.lineY,
          line: { color: '#0ea5e9', width: 2.2 },
          name: 'Líneas de Campo Continuas',
          hoverinfo: 'none',
          showlegend: true
        });
      }

      if (state.showVectors) {
        const quiver = buildQuiverVectors(state.vectorDensity, [-0.95, 0.95], [-0.95, 0.95]);
        data.push({
          type: 'scatter',
          mode: 'lines',
          x: quiver.lineX,
          y: quiver.lineY,
          line: { color: '#38bdf8', width: 1.6 },
          name: 'Campo E = -∇V',
          hoverinfo: 'none',
          showlegend: false
        });
      }

      // Distance line between selected Charge A and B
      const cA = state.charges.find(c => c.id === state.selectedChargeA) || state.charges[0];
      const cB = state.charges.find(c => c.id === state.selectedChargeB) || (state.charges[1] || state.charges[0]);
      if (cA && cB && cA.id !== cB.id) {
        data.push({
          type: 'scatter',
          mode: 'lines',
          x: [cA.x, cB.x],
          y: [cA.y, cB.y],
          line: { color: '#f59e0b', width: 2, dash: 'dot' },
          name: 'Separación r_AB',
          hoverinfo: 'text',
          text: `Separación r_AB`
        });
      }

      // Charge Markers in 2D
      data.push({
        type: 'scatter',
        mode: 'markers+text',
        x: state.charges.map(c => c.x),
        y: state.charges.map(c => c.y),
        text: state.charges.map((c, idx) => `q${idx + 1} (${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC)`),
        textposition: 'top center',
        textfont: { family: 'Inter', size: 11, color: '#ffffff' },
        marker: {
          size: 18,
          color: state.charges.map(c => c.q > 0 ? '#f43f5e' : (c.q < 0 ? '#06b6d4' : '#94a3b8')),
          line: { color: '#ffffff', width: 2.5 }
        },
        name: 'Cargas',
        showlegend: false
      });

      const layout = {
        ...plotlyCommonLayout,
        xaxis: { title: 'x (metros)', range: [-1.05, 1.05], zerolinecolor: '#475569', gridcolor: '#1e293b' },
        yaxis: { title: 'y (metros)', range: [-1.05, 1.05], zerolinecolor: '#475569', gridcolor: '#1e293b', scaleanchor: 'x', scaleratio: 1 }
      };

      if (plotEl._currentMode !== '2d') {
        Plotly.purge(plotEl);
        Plotly.newPlot(plotEl, data, layout, plotlyConfig);
        plotEl._currentMode = '2d';
        attachPlotClickHandler(plotEl);
      } else {
        Plotly.react(plotEl, data, layout, plotlyConfig);
      }
    }
  }

  // --- Module 3 Plot: Gauss's Law 3D Sphere & 2D Cross Section (Robust 3D <-> 2D Purge Transition) ---
  function renderModule3Plot(xArr, yArr, V) {
    const plotEl = elements.plotM3;
    if (!plotEl) return;

    const targetMode = state.m3ViewMode; // '3d' or '2d'
    const data = [];
    const R = state.radius;

    if (targetMode === '3d') {
      // 1. Translucent 3D Gaussian Sphere Surface
      const nU = 26;
      const nV = 16;
      const uArr = linspace(0, 2 * Math.PI, nU);
      const vArr = linspace(0, Math.PI, nV);

      const sphereX = [];
      const sphereY = [];
      const sphereZ = [];

      for (let j = 0; j < nV; j++) {
        const rowX = [];
        const rowY = [];
        const rowZ = [];
        const phi = vArr[j];
        for (let i = 0; i < nU; i++) {
          const theta = uArr[i];
          rowX.push(R * Math.sin(phi) * Math.cos(theta));
          rowY.push(R * Math.sin(phi) * Math.sin(theta));
          rowZ.push(R * Math.cos(phi));
        }
        sphereX.push(rowX);
        sphereY.push(rowY);
        sphereZ.push(rowZ);
      }

      data.push({
        type: 'surface',
        x: sphereX,
        y: sphereY,
        z: sphereZ,
        opacity: 0.22,
        colorscale: [[0, '#10b981'], [1, '#34d399']],
        showscale: false,
        hoverinfo: 'none'
      });

      // 2. Outward Normal Vectors dA along sphere
      const normSample = 14;
      const normLineX = [], normLineY = [], normLineZ = [];
      for (let k = 0; k < normSample; k++) {
        const theta = (k / normSample) * 2 * Math.PI;
        const phi = (k % 2 === 0) ? Math.PI / 3 : (2 * Math.PI / 3);
        const x0 = R * Math.sin(phi) * Math.cos(theta);
        const y0 = R * Math.sin(phi) * Math.sin(theta);
        const z0 = R * Math.cos(phi);

        const len = 0.14;
        const x1 = x0 + (x0 / R) * len;
        const y1 = y0 + (y0 / R) * len;
        const z1 = z0 + (z0 / R) * len;

        normLineX.push(x0, x1, null);
        normLineY.push(y0, y1, null);
        normLineZ.push(z0, z1, null);
      }

      data.push({
        type: 'scatter3d',
        mode: 'lines',
        x: normLineX,
        y: normLineY,
        z: normLineZ,
        line: { color: '#a7f3d0', width: 4 },
        name: 'Vectores Normales dA',
        hoverinfo: 'none'
      });

      // 3. Electric field vectors E at sample surface points
      const efLineX = [], efLineY = [], efLineZ = [];
      for (let k = 0; k < normSample; k++) {
        const theta = (k / normSample) * 2 * Math.PI;
        const x0 = R * Math.cos(theta);
        const y0 = R * Math.sin(theta);
        const z0 = 0;

        const { ex, ey, ez, mag } = evaluateElectricFieldAt(x0, y0, z0);
        if (mag > 1e-4) {
          const len = 0.16;
          const x1 = x0 + (ex / mag) * len;
          const y1 = y0 + (ey / mag) * len;
          const z1 = z0 + (ez / mag) * len;

          efLineX.push(x0, x1, null);
          efLineY.push(y0, y1, null);
          efLineZ.push(z0, z1, null);
        }
      }

      data.push({
        type: 'scatter3d',
        mode: 'lines',
        x: efLineX,
        y: efLineY,
        z: efLineZ,
        line: { color: '#38bdf8', width: 3 },
        name: 'Campo E en superficie',
        hoverinfo: 'none'
      });

      // 3D Continuous Field Streamlines with directional arrows at z = 0
      const streamlines3D = buildFieldStreamlines();
      if (streamlines3D.lineX.length > 0) {
        data.push({
          type: 'scatter3d',
          mode: 'lines',
          x: streamlines3D.lineX,
          y: streamlines3D.lineY,
          z: streamlines3D.lineZ,
          line: { color: '#38bdf8', width: 3 },
          name: 'Líneas de Campo Continuas',
          hoverinfo: 'none',
          showlegend: true
        });
      }

      // 4. All Charges in 3D with Inside/Outside status
      const chargeX = [];
      const chargeY = [];
      const chargeZ = [];
      const chargeText = [];
      const chargeColors = [];

      state.charges.forEach((c, idx) => {
        const dist = Math.sqrt(c.x * c.x + c.y * c.y);
        const isInside = dist <= R;

        chargeX.push(c.x);
        chargeY.push(c.y);
        chargeZ.push(0.0);
        chargeText.push(`q${idx + 1} (${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC)\n[${isInside ? 'DENTRO' : 'FUERA'}]`);
        chargeColors.push(isInside ? (c.q >= 0 ? '#f43f5e' : '#06b6d4') : 'rgba(148, 163, 184, 0.4)');
      });

      data.push({
        type: 'scatter3d',
        mode: 'markers+text',
        x: chargeX,
        y: chargeY,
        z: chargeZ,
        text: chargeText,
        textposition: 'top center',
        textfont: { family: 'Inter', size: 10, color: '#ffffff' },
        marker: {
          size: 9,
          color: chargeColors,
          line: { color: '#ffffff', width: 2 }
        },
        showlegend: false
      });

      const layout = {
        ...plotlyCommonLayout,
        scene: {
          dragmode: 'orbit',
          xaxis: { title: 'x (m)', color: '#94a3b8', gridcolor: '#334155', backgroundcolor: 'rgba(15, 23, 42, 0.6)' },
          yaxis: { title: 'y (m)', color: '#94a3b8', gridcolor: '#334155', backgroundcolor: 'rgba(15, 23, 42, 0.6)' },
          zaxis: { title: 'z (m)', color: '#94a3b8', gridcolor: '#334155', backgroundcolor: 'rgba(15, 23, 42, 0.6)' },
          camera: { eye: { x: 1.5, y: -1.5, z: 1.2 } }
        }
      };

      if (plotEl._currentMode !== '3d') {
        Plotly.purge(plotEl);
        Plotly.newPlot(plotEl, data, layout, plotlyConfig);
        plotEl._currentMode = '3d';
        attachPlotClickHandler(plotEl);
      } else {
        Plotly.react(plotEl, data, layout, plotlyConfig);
      }
    } else {
      // 2D Cross Section View
      // Continuous Field Streamlines with directional arrows in 2D
      const streamlines2D = buildFieldStreamlines();
      if (streamlines2D.lineX.length > 0) {
        data.push({
          type: 'scatter',
          mode: 'lines',
          x: streamlines2D.lineX,
          y: streamlines2D.lineY,
          line: { color: '#0ea5e9', width: 2.2 },
          name: 'Líneas de Campo Continuas',
          hoverinfo: 'none',
          showlegend: true
        });
      }

      const quiver = buildQuiverVectors(17, [-0.95, 0.95], [-0.95, 0.95]);
      data.push({
        type: 'scatter',
        mode: 'lines',
        x: quiver.lineX,
        y: quiver.lineY,
        line: { color: 'rgba(56, 189, 248, 0.45)', width: 1.4 },
        name: 'Líneas de Campo E',
        hoverinfo: 'none',
        showlegend: false
      });

      // Circle Perimeter
      const numPts = 100;
      const circleX = [];
      const circleY = [];
      for (let k = 0; k <= numPts; k++) {
        const theta = (k / numPts) * 2 * Math.PI;
        circleX.push(R * Math.cos(theta));
        circleY.push(R * Math.sin(theta));
      }

      data.push({
        type: 'scatter',
        mode: 'lines',
        x: circleX,
        y: circleY,
        fill: 'toself',
        fillcolor: 'rgba(16, 185, 129, 0.12)',
        line: { color: '#10b981', width: 3, dash: 'dash' },
        name: `Círculo Gaussiano (R = ${R.toFixed(2)} m)`,
        hoverinfo: 'name'
      });

      // Normal vectors dA
      const normalCount = 12;
      const normLineX = [];
      const normLineY = [];
      const normLen = 0.08;
      for (let k = 0; k < normalCount; k++) {
        const theta = (k / normalCount) * 2 * Math.PI;
        const x0 = R * Math.cos(theta);
        const y0 = R * Math.sin(theta);
        const x1 = (R + normLen) * Math.cos(theta);
        const y1 = (R + normLen) * Math.sin(theta);

        normLineX.push(x0, x1, null);
        normLineY.push(y0, y1, null);
      }

      data.push({
        type: 'scatter',
        mode: 'lines',
        x: normLineX,
        y: normLineY,
        line: { color: '#a7f3d0', width: 2 },
        name: 'Vectores Normales dA',
        hoverinfo: 'none'
      });

      // Charge Markers in 2D
      const cX = [], cY = [], cText = [], cColors = [];
      state.charges.forEach((c, idx) => {
        const dist = Math.sqrt(c.x * c.x + c.y * c.y);
        const isInside = dist <= R;

        cX.push(c.x);
        cY.push(c.y);
        cText.push(`q${idx + 1} (${c.q >= 0 ? '+' : ''}${c.q.toFixed(1)} nC)\n[${isInside ? 'DENTRO' : 'FUERA'}]`);
        cColors.push(isInside ? (c.q >= 0 ? '#f43f5e' : '#06b6d4') : 'rgba(148, 163, 184, 0.45)');
      });

      data.push({
        type: 'scatter',
        mode: 'markers+text',
        x: cX,
        y: cY,
        text: cText,
        textposition: 'top center',
        textfont: { family: 'Inter', size: 10, color: '#f8fafc' },
        marker: {
          size: 18,
          color: cColors,
          line: { color: '#ffffff', width: 2.5 }
        },
        name: 'Cargas Eléctricas',
        showlegend: false
      });

      const layout = {
        ...plotlyCommonLayout,
        xaxis: { title: 'x (metros)', range: [-1.05, 1.05], zerolinecolor: '#475569', gridcolor: '#1e293b' },
        yaxis: { title: 'y (metros)', range: [-1.05, 1.05], zerolinecolor: '#475569', gridcolor: '#1e293b', scaleanchor: 'x', scaleratio: 1 }
      };

      if (plotEl._currentMode !== '2d') {
        Plotly.purge(plotEl);
        Plotly.newPlot(plotEl, data, layout, plotlyConfig);
        plotEl._currentMode = '2d';
        attachPlotClickHandler(plotEl);
      } else {
        Plotly.react(plotEl, data, layout, plotlyConfig);
      }
    }
  }

  // --- CSV Export & Import Features ---

  /**
   * Generates CSV string of the grid mesh with complete charge parameters in header
   */
  function generateCSVData() {
    const N = state.resolution;
    const xArr = linspace(-1.0, 1.0, N);
    const yArr = linspace(-1.0, 1.0, N);
    const { V } = computePotentialGrid(xArr, yArr);

    let csv = `# GradVolt Web 3D - Sears & Zemansky Electrostatic Mesh Export\n`;
    csv += `# Parameters: R=${state.radius},N=${N},charges=${JSON.stringify(state.charges)}\n`;
    csv += `x,y,V\n`;

    for (let j = 0; j < N; j++) {
      const y = yArr[j].toFixed(4);
      for (let i = 0; i < N; i++) {
        const x = xArr[i].toFixed(4);
        const v = V[j][i].toFixed(6);
        csv += `${x},${y},${v}\n`;
      }
    }

    return csv;
  }

  function exportCSV() {
    const csvContent = generateCSVData();
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `gradvolt_potencial_${state.charges.length}cargas_${state.resolution}x${state.resolution}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);

    showToast('Malla y parámetros de cargas exportados en formato CSV', 'success');
  }

  function importCSV(file) {
    if (!file) return;

    const reader = new FileReader();
    reader.onload = function (e) {
      try {
        const text = e.target.result;
        const lines = text.split(/\r?\n/).map(l => l.trim()).filter(l => l.length > 0);

        if (lines.length < 2) {
          throw new Error('El archivo CSV está vacío o incompleto.');
        }

        let metaParsed = false;
        for (let i = 0; i < Math.min(6, lines.length); i++) {
          if (lines[i].startsWith('#') && lines[i].includes('charges=')) {
            const paramStr = lines[i].replace('#', '').replace('Parameters:', '').trim();
            const chargesMatch = paramStr.match(/charges=(\[.*?\])/);
            if (chargesMatch && chargesMatch[1]) {
              const loadedCharges = JSON.parse(chargesMatch[1]);
              if (Array.isArray(loadedCharges) && loadedCharges.length > 0) {
                state.charges = loadedCharges.map(c => ({
                  id: c.id,
                  q: parseFloat(c.q),
                  x: parseFloat(c.x),
                  y: parseFloat(c.y),
                  z: 0.00
                }));
                state.nextChargeId = Math.max(...state.charges.map(c => c.id || 0)) + 1;
                metaParsed = true;
              }
            }

            const rMatch = paramStr.match(/R=([0-9.]+)/);
            if (rMatch) state.radius = parseFloat(rMatch[1]);

            const nMatch = paramStr.match(/N=([0-9]+)/);
            if (nMatch) state.resolution = parseInt(nMatch[1], 10);
            break;
          } else if (lines[i].startsWith('#') && lines[i].includes('q1=')) {
            // Backward compatibility
            const parts = lines[i].replace('#', '').replace('Parameters:', '').trim().split(',');
            let lq1 = 5, lx1 = -0.4, ly1 = 0, lq2 = -5, lx2 = 0.4, ly2 = 0;
            parts.forEach(p => {
              const [k, v] = p.split('=').map(s => s.trim());
              const val = parseFloat(v);
              if (!isNaN(val)) {
                if (k === 'q1') lq1 = val;
                if (k === 'x1') lx1 = val;
                if (k === 'y1') ly1 = val;
                if (k === 'q2') lq2 = val;
                if (k === 'x2') lx2 = val;
                if (k === 'y2') ly2 = val;
                if (k === 'R') state.radius = val;
                if (k === 'N') state.resolution = Math.round(val);
              }
            });
            state.charges = [
              { id: 1, q: lq1, x: lx1, y: ly1, z: 0 },
              { id: 2, q: lq2, x: lx2, y: ly2, z: 0 }
            ];
            state.nextChargeId = 3;
            metaParsed = true;
            break;
          }
        }

        renderChargesList();
        updateCoulombDropdowns();
        if (elements.selectResolution) elements.selectResolution.value = String(state.resolution);
        scheduleUpdate();

        showToast(
          metaParsed
            ? `CSV cargado: ${state.charges.length} cargas restauradas exitosamente.`
            : 'Archivo CSV procesado.',
          'success'
        );
      } catch (err) {
        console.error('Error al importar CSV:', err);
        showToast(`Error al importar CSV: ${err.message}`, 'error');
      }
    };
    reader.readAsText(file);
  }

  // --- Update Preview Table in Module 1 ---
  function updatePreviewTable(xArr, yArr, V) {
    if (!elements.csvPreviewBody) return;

    if (elements.tablePointsCount) {
      elements.tablePointsCount.textContent = `(${xArr.length} × ${yArr.length} = ${xArr.length * yArr.length} celdas)`;
    }

    let html = '';
    const sampleLimit = Math.min(15, xArr.length);

    for (let i = 0; i < sampleLimit; i++) {
      const x = xArr[i].toFixed(4);
      const y = yArr[0].toFixed(4);
      const v = V[0][i].toFixed(4);
      html += `
        <tr class="hover:bg-slate-800/40">
          <td class="px-3 py-1.5 text-slate-500 font-mono">${i + 1}</td>
          <td class="px-3 py-1.5">${x}</td>
          <td class="px-3 py-1.5">${y}</td>
          <td class="px-3 py-1.5 font-semibold ${v >= 0 ? 'text-rose-400' : 'text-cyan-400'}">${v >= 0 ? '+' : ''}${v}</td>
        </tr>
      `;
    }

    elements.csvPreviewBody.innerHTML = html;
  }

  // --- Scenarios ("Escenarios Rápidos para la Feria") ---
  function applyScenario(scenarioId) {
    switch (scenarioId) {
      case 1:
        // Reto 1: Ley Inverso Cuadrado (Carga puntual única en el origen)
        state.charges = [
          { id: 1, q: 8.0, x: 0.00, y: 0.00, z: 0.00 }
        ];
        state.nextChargeId = 2;
        state.radius = 0.40;
        setActiveTab('module-1');
        showToast('Reto 1: Carga puntual única. Analice la caída radial 1/r (V) y 1/r² (E).', 'info');
        break;

      case 2:
        // Reto 2: Paradoja Esfera de Gauss (1 carga dentro, 1 fuera)
        state.charges = [
          { id: 1, q: 5.0, x: 0.00, y: 0.00, z: 0.00 },
          { id: 2, q: -5.0, x: 0.70, y: 0.70, z: 0.00 }
        ];
        state.nextChargeId = 3;
        state.radius = 0.50;
        setActiveTab('module-3');
        showToast('Reto 2: ¡El flujo Φ_E no cambia al variar R mientras q₂ siga fuera!', 'info');
        break;

      case 3:
        // Reto 3: Dipolo Eléctrico Simétrico (V=0, E ≠ 0 en el plano medio)
        state.charges = [
          { id: 1, q: 6.0, x: -0.40, y: 0.00, z: 0.00 },
          { id: 2, q: -6.0, x: 0.40, y: 0.00, z: 0.00 }
        ];
        state.nextChargeId = 3;
        state.radius = 0.60;
        setActiveTab('module-2');
        showToast('Reto 3: En el eje de simetría x=0, V=0 V pero el campo E es perpendicular y máximo.', 'info');
        break;
    }

    state.placingCharge = false;
    updateAddButtonUI();
    renderChargesList();
    updateCoulombDropdowns();
    scheduleUpdate();
  }

  // --- Tab Management ---
  function setActiveTab(tabId) {
    state.activeTab = tabId;

    document.querySelectorAll('.tab-btn').forEach(btn => {
      btn.classList.toggle('active', btn.getAttribute('data-tab') === tabId);
    });

    document.querySelectorAll('.module-tab-content').forEach(section => {
      section.classList.toggle('hidden', section.id !== tabId);
    });

    setTimeout(() => {
      if (tabId === 'module-1' && elements.plotM1) Plotly.Plots.resize(elements.plotM1);
      if (tabId === 'module-2' && elements.plotM2) Plotly.Plots.resize(elements.plotM2);
      if (tabId === 'module-3' && elements.plotM3) Plotly.Plots.resize(elements.plotM3);
      scheduleUpdate();
    }, 60);
  }

  // --- Toast Notification System ---
  function showToast(message, type = 'info') {
    if (!elements.toastContainer) return;

    const toast = document.createElement('div');
    toast.className = `toast px-4 py-2.5 rounded-xl text-xs font-semibold shadow-xl border flex items-center gap-2 pointer-events-auto transition-all ${
      type === 'success'
        ? 'bg-emerald-950/90 text-emerald-200 border-emerald-500/50'
        : type === 'error'
        ? 'bg-rose-950/90 text-rose-200 border-rose-500/50'
        : 'bg-slate-900/90 text-sky-200 border-sky-500/50'
    }`;

    const icon = type === 'success' ? 'fa-circle-check text-emerald-400' : (type === 'error' ? 'fa-triangle-exclamation text-rose-400' : 'fa-circle-info text-sky-400');
    toast.innerHTML = `<i class="fa-solid ${icon}"></i><span>${message}</span>`;

    elements.toastContainer.appendChild(toast);

    setTimeout(() => {
      toast.style.opacity = '0';
      toast.style.transform = 'translateY(10px)';
      setTimeout(() => toast.remove(), 300);
    }, 3800);
  }

  // --- View Mode Button Helpers ---
  function updateM1Buttons() {
    if (!elements.btnView3D || !elements.btnViewContour) return;
    if (state.m1ViewMode === '3d') {
      elements.btnView3D.className = 'px-3 py-1 rounded-lg text-xs font-medium bg-brand-600 text-white shadow transition-all flex items-center gap-1.5';
      elements.btnViewContour.className = 'px-3 py-1 rounded-lg text-xs font-medium text-slate-400 hover:text-white transition-all flex items-center gap-1.5';
    } else {
      elements.btnViewContour.className = 'px-3 py-1 rounded-lg text-xs font-medium bg-brand-600 text-white shadow transition-all flex items-center gap-1.5';
      elements.btnView3D.className = 'px-3 py-1 rounded-lg text-xs font-medium text-slate-400 hover:text-white transition-all flex items-center gap-1.5';
    }
  }

  function updateM2Buttons() {
    if (!elements.btnM2View3D || !elements.btnM2View2D) return;
    if (state.m2ViewMode === '3d') {
      elements.btnM2View3D.className = 'px-3 py-1 rounded-lg text-xs font-medium bg-brand-600 text-white shadow transition-all flex items-center gap-1.5';
      elements.btnM2View2D.className = 'px-3 py-1 rounded-lg text-xs font-medium text-slate-400 hover:text-white transition-all flex items-center gap-1.5';
    } else {
      elements.btnM2View2D.className = 'px-3 py-1 rounded-lg text-xs font-medium bg-brand-600 text-white shadow transition-all flex items-center gap-1.5';
      elements.btnM2View3D.className = 'px-3 py-1 rounded-lg text-xs font-medium text-slate-400 hover:text-white transition-all flex items-center gap-1.5';
    }
  }

  function updateM3Buttons() {
    if (!elements.btnM3View3D || !elements.btnM3View2D) return;
    if (state.m3ViewMode === '3d') {
      elements.btnM3View3D.className = 'px-3 py-1 rounded-lg text-xs font-medium bg-emerald-600 text-white shadow transition-all flex items-center gap-1.5';
      elements.btnM3View2D.className = 'px-3 py-1 rounded-lg text-xs font-medium text-slate-400 hover:text-white transition-all flex items-center gap-1.5';
    } else {
      elements.btnM3View2D.className = 'px-3 py-1 rounded-lg text-xs font-medium bg-emerald-600 text-white shadow transition-all flex items-center gap-1.5';
      elements.btnM3View3D.className = 'px-3 py-1 rounded-lg text-xs font-medium text-slate-400 hover:text-white transition-all flex items-center gap-1.5';
    }
  }

  // --- Event Listeners Setup (Clean, No Uninitialized Plotly .on Calls) ---
  function setupEventListeners() {
    // 1. Direct Add Charge Button
    if (elements.btnAddCharge) {
      elements.btnAddCharge.addEventListener('click', handleAddChargeClick);
    }

    // 2. Gauss Radius Slider
    if (elements.sliderRadius) {
      elements.sliderRadius.addEventListener('input', e => {
        state.radius = parseFloat(e.target.value);
        scheduleUpdate();
      });
    }

    // 3. Resolution Selector
    if (elements.selectResolution) {
      elements.selectResolution.addEventListener('change', e => {
        state.resolution = parseInt(e.target.value, 10);
        scheduleUpdate();
      });
    }

    // 4. Reset Defaults
    if (elements.btnResetDefaults) {
      elements.btnResetDefaults.addEventListener('click', () => {
        state.charges = [
          { id: 1, q: 5.0, x: -0.40, y: 0.00, z: 0.00 },
          { id: 2, q: -5.0, x: 0.40, y: 0.00, z: 0.00 }
        ];
        state.nextChargeId = 3;
        state.radius = 0.50;
        state.resolution = 41;
        state.selectedChargeA = 1;
        state.selectedChargeB = 2;
        state.placingCharge = false;
        state.importedFromCsv = false;
        if (elements.selectResolution) elements.selectResolution.value = '41';
        renderChargesList();
        updateCoulombDropdowns();
        updateAddButtonUI();
        scheduleUpdate();
        showToast('Parámetros restablecidos a los valores por defecto.', 'info');
      });
    }

    // 5. Module 1 View Mode toggles with Plotly purge
    if (elements.btnView3D) {
      elements.btnView3D.addEventListener('click', () => {
        if (state.m1ViewMode !== '3d') {
          state.m1ViewMode = '3d';
          updateM1Buttons();
          if (elements.plotM1) {
            Plotly.purge(elements.plotM1);
            elements.plotM1._currentMode = null;
          }
          scheduleUpdate();
        }
      });
    }

    if (elements.btnViewContour) {
      elements.btnViewContour.addEventListener('click', () => {
        if (state.m1ViewMode !== 'contour') {
          state.m1ViewMode = 'contour';
          updateM1Buttons();
          if (elements.plotM1) {
            Plotly.purge(elements.plotM1);
            elements.plotM1._currentMode = null;
          }
          scheduleUpdate();
        }
      });
    }

    if (elements.selectColorscale) {
      elements.selectColorscale.addEventListener('change', e => {
        state.colorscale = e.target.value;
        scheduleUpdate();
      });
    }

    // Module 1 Table Toggle
    if (elements.toggleTableBtn && elements.tableContainer) {
      elements.toggleTableBtn.addEventListener('click', () => {
        const isHidden = elements.tableContainer.classList.contains('hidden');
        elements.tableContainer.classList.toggle('hidden');
        if (elements.tableChevron) {
          elements.tableChevron.style.transform = isHidden ? 'rotate(180deg)' : 'rotate(0deg)';
        }
      });
    }

    // 6. Module 2 View Mode toggles with Plotly purge
    if (elements.btnM2View3D) {
      elements.btnM2View3D.addEventListener('click', () => {
        if (state.m2ViewMode !== '3d') {
          state.m2ViewMode = '3d';
          updateM2Buttons();
          if (elements.plotM2) {
            Plotly.purge(elements.plotM2);
            elements.plotM2._currentMode = null;
          }
          scheduleUpdate();
        }
      });
    }

    if (elements.btnM2View2D) {
      elements.btnM2View2D.addEventListener('click', () => {
        if (state.m2ViewMode !== '2d') {
          state.m2ViewMode = '2d';
          updateM2Buttons();
          if (elements.plotM2) {
            Plotly.purge(elements.plotM2);
            elements.plotM2._currentMode = null;
          }
          scheduleUpdate();
        }
      });
    }

    // Module 2 Coulomb Pair Selectors
    if (elements.selectCoulombA) {
      elements.selectCoulombA.addEventListener('change', e => {
        state.selectedChargeA = parseInt(e.target.value, 10);
        scheduleUpdate();
      });
    }

    if (elements.selectCoulombB) {
      elements.selectCoulombB.addEventListener('change', e => {
        state.selectedChargeB = parseInt(e.target.value, 10);
        scheduleUpdate();
      });
    }

    if (elements.checkShowEquipotential) {
      elements.checkShowEquipotential.addEventListener('change', e => {
        state.showEquipotential = e.target.checked;
        scheduleUpdate();
      });
    }

    if (elements.checkShowVectors) {
      elements.checkShowVectors.addEventListener('change', e => {
        state.showVectors = e.target.checked;
        scheduleUpdate();
      });
    }

    if (elements.selectVectorDensity) {
      elements.selectVectorDensity.addEventListener('change', e => {
        state.vectorDensity = parseInt(e.target.value, 10);
        scheduleUpdate();
      });
    }

    // 7. Module 3 View Mode toggles with Plotly purge
    if (elements.btnM3View3D) {
      elements.btnM3View3D.addEventListener('click', () => {
        if (state.m3ViewMode !== '3d') {
          state.m3ViewMode = '3d';
          updateM3Buttons();
          if (elements.plotM3) {
            Plotly.purge(elements.plotM3);
            elements.plotM3._currentMode = null;
          }
          scheduleUpdate();
        }
      });
    }

    if (elements.btnM3View2D) {
      elements.btnM3View2D.addEventListener('click', () => {
        if (state.m3ViewMode !== '2d') {
          state.m3ViewMode = '2d';
          updateM3Buttons();
          if (elements.plotM3) {
            Plotly.purge(elements.plotM3);
            elements.plotM3._currentMode = null;
          }
          scheduleUpdate();
        }
      });
    }

    // 8. Navigation Tab buttons
    document.querySelectorAll('.tab-btn').forEach(btn => {
      btn.addEventListener('click', () => {
        const tab = btn.getAttribute('data-tab');
        setActiveTab(tab);
      });
    });

    // 9. Quick Fair Scenarios (Retos)
    if (elements.btnScenario1) elements.btnScenario1.addEventListener('click', () => applyScenario(1));
    if (elements.btnScenario2) elements.btnScenario2.addEventListener('click', () => applyScenario(2));
    if (elements.btnScenario3) elements.btnScenario3.addEventListener('click', () => applyScenario(3));

    // 10. CSV Export & Import
    if (elements.btnExportCsv) elements.btnExportCsv.addEventListener('click', exportCSV);
    if (elements.btnExportCsvSub) elements.btnExportCsvSub.addEventListener('click', exportCSV);

    if (elements.inputImportCsv) {
      elements.inputImportCsv.addEventListener('change', e => {
        const file = e.target.files[0];
        if (file) {
          importCSV(file);
          e.target.value = '';
        }
      });
    }

    // 11. Theory Modal
    if (elements.btnOpenTheory && elements.modalTheory) {
      elements.btnOpenTheory.addEventListener('click', () => {
        elements.modalTheory.classList.remove('hidden');
        renderMathInElement();
      });
    }
    const closeModal = () => {
      if (elements.modalTheory) elements.modalTheory.classList.add('hidden');
    };
    if (elements.btnCloseTheory) elements.btnCloseTheory.addEventListener('click', closeModal);
    if (elements.btnCloseTheory2) elements.btnCloseTheory2.addEventListener('click', closeModal);
    if (elements.modalTheory) {
      elements.modalTheory.addEventListener('click', e => {
        if (e.target === elements.modalTheory) closeModal();
      });
    }

    // 12. Window Resize
    window.addEventListener('resize', () => {
      if (state.activeTab === 'module-1' && elements.plotM1) Plotly.Plots.resize(elements.plotM1);
      if (state.activeTab === 'module-2' && elements.plotM2) Plotly.Plots.resize(elements.plotM2);
      if (state.activeTab === 'module-3' && elements.plotM3) Plotly.Plots.resize(elements.plotM3);
    });
  }

  // --- KaTeX Auto-Render Helper ---
  function renderMathInElement() {
    if (window.renderMathInElement) {
      try {
        window.renderMathInElement(document.body, {
          delimiters: [
            { left: '$$', right: '$$', display: true },
            { left: '$', right: '$', display: false }
          ],
          throwOnError: false
        });
      } catch (err) {
        console.warn('KaTeX rendering notice:', err);
      }
    }
  }

  // --- Initialization ---
  function init() {
    renderChargesList();
    setupEventListeners();
    scheduleUpdate();

    setTimeout(() => {
      renderMathInElement();
    }, 400);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }

})();
