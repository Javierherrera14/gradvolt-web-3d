/**
 * ============================================================================
 * GradVolt Web 3D - Computational Physics Simulation Platform
 * Based on University Physics with Modern Physics by Sears & Zemansky
 * Chapters 21 (Coulomb & Electric Field), 22 (Gauss's Law), 23 (Potential)
 * ============================================================================
 */

(function () {
  'use strict';

  // --- Physical Constants ---
  const K_COULOMB = 8.9875517923e9; // N·m²/C² (k = 1 / (4πε₀))
  const EPSILON_0 = 8.8541878176e-12; // C²/(N·m²)
  const EPS_SMOOTH = 0.06; // Softening parameter (meters) to prevent 1/0 singularities

  // --- Simulation State ---
  const state = {
    q1: 5.0, // nC
    x1: -0.40, // m
    y1: 0.00, // m
    q2: -5.0, // nC
    x2: 0.40, // m
    y2: 0.00, // m
    radius: 0.50, // m (Gaussian surface radius)
    resolution: 41, // N x N mesh grid
    activeTab: 'module-1',
    m1ViewMode: '3d', // '3d' or 'contour'
    colorscale: 'Electric',
    showEquipotential: true,
    showVectors: true,
    vectorDensity: 19,
    importedFromCsv: false,
    csvMesh: null // Custom mesh if imported
  };

  // --- Cached DOM Elements ---
  const elements = {
    // Sliders & Values
    sliderQ1: document.getElementById('slider-q1'),
    q1Val: document.getElementById('q1-val'),
    q1Badge: document.getElementById('q1-badge'),
    sliderX1: document.getElementById('slider-x1'),
    x1Val: document.getElementById('x1-val'),
    sliderY1: document.getElementById('slider-y1'),
    y1Val: document.getElementById('y1-val'),

    sliderQ2: document.getElementById('slider-q2'),
    q2Val: document.getElementById('q2-val'),
    q2Badge: document.getElementById('q2-badge'),
    sliderX2: document.getElementById('slider-x2'),
    x2Val: document.getElementById('x2-val'),
    sliderY2: document.getElementById('slider-y2'),
    y2Val: document.getElementById('y2-val'),

    sliderRadius: document.getElementById('slider-radius'),
    radiusVal: document.getElementById('radius-val'),
    gaussRBadge: document.getElementById('gauss-r-badge'),
    selectResolution: document.getElementById('select-resolution'),

    // Module 1 controls & metrics
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

    // Module 2 controls & metrics
    coulombDist: document.getElementById('coulomb-dist'),
    coulombForce: document.getElementById('coulomb-force'),
    coulombNature: document.getElementById('coulomb-nature'),
    coulombNatureIcon: document.getElementById('coulomb-nature-icon'),
    checkShowEquipotential: document.getElementById('check-show-equipotential'),
    checkShowVectors: document.getElementById('check-show-vectors'),
    selectVectorDensity: document.getElementById('select-vector-density'),
    plotM2: document.getElementById('plot-module-2'),

    // Module 3 controls & metrics
    gaussQenc: document.getElementById('gauss-qenc'),
    gaussQencSi: document.getElementById('gauss-qenc-si'),
    gaussFlux: document.getElementById('gauss-flux'),
    gaussArea: document.getElementById('gauss-area'),
    gaussRadiusDisp: document.getElementById('gauss-radius-disp'),
    gaussEavg: document.getElementById('gauss-eavg'),
    badgeCharge1Status: document.getElementById('badge-charge-1-status'),
    badgeCharge2Status: document.getElementById('badge-charge-2-status'),
    gaussInvarianceMsg: document.getElementById('gauss-invariance-msg'),
    plotM3: document.getElementById('plot-module-3'),

    // CSV & Scenarios
    btnExportCsv: document.getElementById('btn-export-csv'),
    btnExportCsvSub: document.getElementById('btn-export-csv-sub'),
    inputImportCsv: document.getElementById('input-import-csv'),
    btnResetDefaults: document.getElementById('btn-reset-defaults'),
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
   * Generates a 1D uniform grid from min to max with N points
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
   * Evaluates electric potential at a point (x, y) given charges q1, q2
   * V(x, y) = k * q1 / r1 + k * q2 / r2
   */
  function evaluatePotentialAt(x, y, q1, x1, y1, q2, x2, y2) {
    const q1_SI = q1 * 1e-9;
    const q2_SI = q2 * 1e-9;

    const r1 = Math.sqrt((x - x1) ** 2 + (y - y1) ** 2 + EPS_SMOOTH ** 2);
    const r2 = Math.sqrt((x - x2) ** 2 + (y - y2) ** 2 + EPS_SMOOTH ** 2);

    return K_COULOMB * (q1_SI / r1 + q2_SI / r2);
  }

  /**
   * Generates scalar potential matrix V[j][i] over 2D mesh
   */
  function computePotentialGrid(xArr, yArr, q1, x1, y1, q2, x2, y2) {
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
        const val = evaluatePotentialAt(x, y, q1, x1, y1, q2, x2, y2);
        row[i] = val;
        if (val < vMin) vMin = val;
        if (val > vMax) vMax = val;
      }
      V[j] = row;
    }

    return { V, vMin, vMax };
  }

  /**
   * Computes numerical gradient E = -∇V using Central Finite Differences
   * Ex = - (V[j][i+1] - V[j][i-1]) / (2 * Δx)
   * Ey = - (V[j+1][i] - V[j-1][i]) / (2 * Δy)
   */
  function computeNumericalGradient(V, xArr, yArr) {
    const Ny = yArr.length;
    const Nx = xArr.length;
    const dx = xArr[1] - xArr[0];
    const dy = yArr[1] - yArr[0];

    const Ex = new Array(Ny);
    const Ey = new Array(Ny);
    const Emag = new Array(Ny);

    for (let j = 0; j < Ny; j++) {
      Ex[j] = new Float64Array(Nx);
      Ey[j] = new Float64Array(Nx);
      Emag[j] = new Float64Array(Nx);

      for (let i = 0; i < Nx; i++) {
        // Ex: Central difference, forward/backward at borders
        let dVdx;
        if (i === 0) {
          dVdx = (V[j][1] - V[j][0]) / dx;
        } else if (i === Nx - 1) {
          dVdx = (V[j][Nx - 1] - V[j][Nx - 2]) / dx;
        } else {
          dVdx = (V[j][i + 1] - V[j][i - 1]) / (2 * dx);
        }

        // Ey: Central difference, forward/backward at borders
        let dVdy;
        if (j === 0) {
          dVdy = (V[1][i] - V[0][i]) / dy;
        } else if (j === Ny - 1) {
          dVdy = (V[Ny - 1][i] - V[Ny - 2][i]) / dy;
        } else {
          dVdy = (V[j + 1][i] - V[j - 1][i]) / (2 * dy);
        }

        const ex = -dVdx;
        const ey = -dVdy;
        Ex[j][i] = ex;
        Ey[j][i] = ey;
        Emag[j][i] = Math.sqrt(ex * ex + ey * ey);
      }
    }

    return { Ex, Ey, Emag };
  }

  /**
   * Formats numbers into scientific or human notation
   */
  function formatSci(val, unit = '', decimals = 2) {
    if (Math.abs(val) === 0) return `0.00 ${unit}`.trim();
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
      '-': '⁻',
      '0': '⁰',
      '1': '¹',
      '2': '²',
      '3': '³',
      '4': '⁴',
      '5': '⁵',
      '6': '⁶',
      '7': '⁷',
      '8': '⁸',
      '9': '⁹'
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
      const potResult = computePotentialGrid(xArr, yArr, state.q1, state.x1, state.y1, state.q2, state.x2, state.y2);
      V = potResult.V;
      vMin = potResult.vMin;
      vMax = potResult.vMax;
    }

    // 2. Metrics for Module 1
    const vOrigin = evaluatePotentialAt(0, 0, state.q1, state.x1, state.y1, state.q2, state.x2, state.y2);
    elements.m1Vmax.textContent = `${vMax >= 0 ? '+' : ''}${formatSci(vMax, 'V')}`;
    elements.m1Vmin.textContent = `${formatSci(vMin, 'V')}`;
    elements.m1Vorigin.textContent = `${formatSci(vOrigin, 'V')}`;

    // 3. Coulomb Calculation (Module 2)
    const r12 = Math.sqrt((state.x1 - state.x2) ** 2 + (state.y1 - state.y2) ** 2);
    const q1_SI = state.q1 * 1e-9;
    const q2_SI = state.q2 * 1e-9;
    const forceCoulomb = (K_COULOMB * Math.abs(q1_SI * q2_SI)) / Math.max(r12 * r12, 1e-6);

    elements.coulombDist.textContent = `${r12.toFixed(3)} m`;
    elements.coulombForce.textContent = formatSci(forceCoulomb, 'N');

    if (state.q1 === 0 || state.q2 === 0) {
      elements.coulombNature.textContent = 'Nula (Una o ambas cargas son 0)';
      elements.coulombNature.className = 'text-xs font-semibold text-slate-400';
      elements.coulombNatureIcon.className = 'w-10 h-10 rounded-lg bg-slate-800 border border-slate-700 flex items-center justify-center text-slate-400 text-lg';
    } else if (state.q1 * state.q2 > 0) {
      elements.coulombNature.textContent = 'Repulsiva (Signos iguales)';
      elements.coulombNature.className = 'text-xs font-semibold text-rose-400';
      elements.coulombNatureIcon.className = 'w-10 h-10 rounded-lg bg-rose-500/10 border border-rose-500/20 flex items-center justify-center text-rose-400 text-lg';
    } else {
      elements.coulombNature.textContent = 'Atractiva (Signos opuestos)';
      elements.coulombNature.className = 'text-xs font-semibold text-cyan-400';
      elements.coulombNatureIcon.className = 'w-10 h-10 rounded-lg bg-cyan-500/10 border border-cyan-500/20 flex items-center justify-center text-cyan-400 text-lg';
    }

    // 4. Gauss Calculation (Module 3)
    const d1 = Math.sqrt(state.x1 ** 2 + state.y1 ** 2);
    const d2 = Math.sqrt(state.x2 ** 2 + state.y2 ** 2);
    const q1Inside = d1 <= state.radius;
    const q2Inside = d2 <= state.radius;

    let qEnc = 0;
    if (q1Inside) qEnc += state.q1;
    if (q2Inside) qEnc += state.q2;

    const qEncSI = qEnc * 1e-9;
    const fluxTheoretical = qEncSI / EPSILON_0;
    const areaSphere = 4 * Math.PI * (state.radius ** 2);
    const eAvg = Math.abs(fluxTheoretical) / areaSphere;

    elements.gaussQenc.textContent = `${qEnc >= 0 ? '+' : ''}${qEnc.toFixed(2)} nC`;
    elements.gaussQencSi.textContent = formatSci(qEncSI, 'C');
    elements.gaussFlux.textContent = formatSci(fluxTheoretical, 'N·m²/C');
    elements.gaussArea.textContent = `${areaSphere.toFixed(3)} m²`;
    elements.gaussRadiusDisp.textContent = `R = ${state.radius.toFixed(2)} m`;
    elements.gaussEavg.textContent = formatSci(eAvg, 'N/C');

    // Charge Inside/Outside badges
    elements.badgeCharge1Status.className = `px-2.5 py-1 rounded-full text-[11px] font-semibold flex items-center gap-1.5 ${
      q1Inside 
        ? 'bg-rose-500/20 text-rose-300 border border-rose-500/40 shadow-[0_0_8px_rgba(244,63,94,0.3)]' 
        : 'bg-slate-800 text-slate-500 border border-slate-700'
    }`;
    elements.badgeCharge1Status.innerHTML = `<i class="fa-solid fa-circle text-[8px]"></i><span>q₁: ${q1Inside ? 'DENTRO de la superficie' : 'FUERA de la superficie'}</span>`;

    elements.badgeCharge2Status.className = `px-2.5 py-1 rounded-full text-[11px] font-semibold flex items-center gap-1.5 ${
      q2Inside 
        ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/40 shadow-[0_0_8px_rgba(6,182,212,0.3)]' 
        : 'bg-slate-800 text-slate-500 border border-slate-700'
    }`;
    elements.badgeCharge2Status.innerHTML = `<i class="fa-solid fa-circle text-[8px]"></i><span>q₂: ${q2Inside ? 'DENTRO de la superficie' : 'FUERA de la superficie'}</span>`;

    // Gauss message banner dynamic styling
    if (Math.abs(qEnc) < 1e-6) {
      elements.gaussInvarianceMsg.innerHTML = 'Carga neta interior <strong class="text-amber-300">Q_enc = 0</strong>: El flujo eléctrico total entrante equivale con exactitud al flujo saliente (<strong class="text-white">Φ_E = 0</strong>).';
    } else {
      elements.gaussInvarianceMsg.innerHTML = `Flujo constante en <strong class="text-emerald-300">${formatSci(fluxTheoretical, 'N·m²/C')}</strong>. Al mover el slider de radio R, el flujo permanece rigurosamente idéntico mientras no cruce ninguna carga.`;
    }

    // 5. Update Active Module Visualizations
    if (state.activeTab === 'module-1') {
      renderModule1Plot(xArr, yArr, V);
    } else if (state.activeTab === 'module-2') {
      const gradResult = computeNumericalGradient(V, xArr, yArr);
      renderModule2Plot(xArr, yArr, V, gradResult.Ex, gradResult.Ey);
    } else if (state.activeTab === 'module-3') {
      const gradResult = computeNumericalGradient(V, xArr, yArr);
      renderModule3Plot(xArr, yArr, V, gradResult.Ex, gradResult.Ey);
    }

    // 6. Update CSV Preview Table
    updatePreviewTable(xArr, yArr, V);
  }

  /**
   * Syncs HTML control values and badges with state
   */
  function updateUIElements() {
    elements.sliderQ1.value = state.q1;
    elements.q1Val.textContent = `${state.q1 >= 0 ? '+' : ''}${state.q1.toFixed(1)} nC`;
    elements.q1Badge.textContent = `${state.q1 >= 0 ? '+' : ''}${state.q1.toFixed(1)} nC`;
    elements.sliderX1.value = state.x1;
    elements.x1Val.textContent = `${state.x1.toFixed(2)} m`;
    elements.sliderY1.value = state.y1;
    elements.y1Val.textContent = `${state.y1.toFixed(2)} m`;

    elements.sliderQ2.value = state.q2;
    elements.q2Val.textContent = `${state.q2 >= 0 ? '+' : ''}${state.q2.toFixed(1)} nC`;
    elements.q2Badge.textContent = `${state.q2 >= 0 ? '+' : ''}${state.q2.toFixed(1)} nC`;
    elements.sliderX2.value = state.x2;
    elements.x2Val.textContent = `${state.x2.toFixed(2)} m`;
    elements.sliderY2.value = state.y2;
    elements.y2Val.textContent = `${state.y2.toFixed(2)} m`;

    elements.sliderRadius.value = state.radius;
    elements.radiusVal.textContent = `${state.radius.toFixed(2)} m`;
    elements.gaussRBadge.textContent = `R = ${state.radius.toFixed(2)} m`;

    elements.selectColorscale.value = state.colorscale;
  }

  // --- Module 1 Plot: 3D Surface / 2D Contours ---
  function renderModule1Plot(xArr, yArr, V) {
    if (!elements.plotM1) return;

    const data = [];

    if (state.m1ViewMode === '3d') {
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

      const layout = {
        ...plotlyCommonLayout,
        scene: {
          xaxis: {
            title: 'x (m)',
            color: '#94a3b8',
            gridcolor: '#334155',
            showbackground: true,
            backgroundcolor: 'rgba(15, 23, 42, 0.6)'
          },
          yaxis: {
            title: 'y (m)',
            color: '#94a3b8',
            gridcolor: '#334155',
            showbackground: true,
            backgroundcolor: 'rgba(15, 23, 42, 0.6)'
          },
          zaxis: {
            title: 'Potencial V (V)',
            color: '#94a3b8',
            gridcolor: '#334155',
            showbackground: true,
            backgroundcolor: 'rgba(15, 23, 42, 0.6)'
          },
          camera: {
            eye: { x: 1.45, y: -1.45, z: 1.25 }
          }
        }
      };

      Plotly.react(elements.plotM1, data, layout, plotlyConfig);
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
        line: {
          width: 0.8,
          color: 'rgba(255,255,255,0.4)'
        },
        colorbar: {
          title: { text: 'V (Voltios)', font: { color: '#e2e8f0', size: 12 } },
          tickfont: { color: '#94a3b8', size: 10 },
          thickness: 16
        }
      });

      // Overlay charge positions
      data.push({
        type: 'scatter',
        mode: 'markers+text',
        x: [state.x1, state.x2],
        y: [state.y1, state.y2],
        text: [`q₁ (${state.q1 >= 0 ? '+' : ''}${state.q1} nC)`, `q₂ (${state.q2 >= 0 ? '+' : ''}${state.q2} nC)`],
        textposition: 'top center',
        textfont: { family: 'Inter', size: 11, color: '#f8fafc' },
        marker: {
          size: 16,
          color: [state.q1 >= 0 ? '#f43f5e' : '#06b6d4', state.q2 >= 0 ? '#f43f5e' : '#06b6d4'],
          line: { color: '#ffffff', width: 2 }
        },
        showlegend: false
      });

      const layout = {
        ...plotlyCommonLayout,
        xaxis: {
          title: 'x (metros)',
          range: [-1.05, 1.05],
          zerolinecolor: '#475569',
          gridcolor: '#1e293b'
        },
        yaxis: {
          title: 'y (metros)',
          range: [-1.05, 1.05],
          zerolinecolor: '#475569',
          gridcolor: '#1e293b',
          scaleanchor: 'x',
          scaleratio: 1
        }
      };

      Plotly.react(elements.plotM1, data, layout, plotlyConfig);
    }
  }

  // --- Helper to build vector quiver lines ---
  function buildQuiverVectors(density, xRange, yRange) {
    const qX = linspace(xRange[0], xRange[1], density);
    const qY = linspace(yRange[0], yRange[1], density);
    const arrowLen = 1.4 / density; // Arrow visual scale
    const headLen = arrowLen * 0.35;
    const headAngle = 0.45; // radians

    const lineX = [];
    const lineY = [];

    for (let j = 0; j < density; j++) {
      const y0 = qY[j];
      for (let i = 0; i < density; i++) {
        const x0 = qX[i];

        // Numerical gradient at (x0, y0) using small delta
        const delta = 0.005;
        const vRight = evaluatePotentialAt(x0 + delta, y0, state.q1, state.x1, state.y1, state.q2, state.x2, state.y2);
        const vLeft = evaluatePotentialAt(x0 - delta, y0, state.q1, state.x1, state.y1, state.q2, state.x2, state.y2);
        const vTop = evaluatePotentialAt(x0, y0 + delta, state.q1, state.x1, state.y1, state.q2, state.x2, state.y2);
        const vBottom = evaluatePotentialAt(x0, y0 - delta, state.q1, state.x1, state.y1, state.q2, state.x2, state.y2);

        const ex = -(vRight - vLeft) / (2 * delta);
        const ey = -(vTop - vBottom) / (2 * delta);
        const mag = Math.sqrt(ex * ex + ey * ey);

        if (mag < 1e-4) continue;

        // Normalized unit vector
        const ux = ex / mag;
        const uy = ey / mag;

        // Arrow endpoints (centered at grid point)
        const startX = x0 - (ux * arrowLen * 0.4);
        const startY = y0 - (uy * arrowLen * 0.4);
        const endX = x0 + (ux * arrowLen * 0.6);
        const endY = y0 + (uy * arrowLen * 0.6);

        // Arrowhead barb 1
        const b1X = endX - headLen * (ux * Math.cos(headAngle) - uy * Math.sin(headAngle));
        const b1Y = endY - headLen * (uy * Math.cos(headAngle) + ux * Math.sin(headAngle));

        // Arrowhead barb 2
        const b2X = endX - headLen * (ux * Math.cos(-headAngle) - uy * Math.sin(-headAngle));
        const b2Y = endY - headLen * (uy * Math.cos(-headAngle) + ux * Math.sin(-headAngle));

        // Append line segments separated by null
        lineX.push(startX, endX, b1X, endX, b2X, null);
        lineY.push(startY, endY, b1Y, endY, b2Y, null);
      }
    }

    return { lineX, lineY };
  }

  // --- Module 2 Plot: E-field via -∇V + Equipotentials ---
  function renderModule2Plot(xArr, yArr, V) {
    if (!elements.plotM2) return;

    const data = [];

    // 1. Background Equipotential Contours
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
        colorbar: {
          title: { text: 'V (V)', font: { color: '#e2e8f0', size: 11 } },
          thickness: 14,
          len: 0.8
        },
        hoverinfo: 'x+y+z'
      });
    }

    // 2. Vector field (Quiver arrows)
    if (state.showVectors) {
      const quiver = buildQuiverVectors(state.vectorDensity, [-0.95, 0.95], [-0.95, 0.95]);
      data.push({
        type: 'scatter',
        mode: 'lines',
        x: quiver.lineX,
        y: quiver.lineY,
        line: {
          color: '#38bdf8',
          width: 1.6
        },
        name: 'Campo E = -∇V',
        hoverinfo: 'none',
        showlegend: false
      });
    }

    // 3. Distance line between charges
    data.push({
      type: 'scatter',
      mode: 'lines',
      x: [state.x1, state.x2],
      y: [state.y1, state.y2],
      line: {
        color: '#f59e0b',
        width: 2,
        dash: 'dot'
      },
      name: 'Separación r₁₂',
      hoverinfo: 'text',
      text: `Distancia r₁₂: ${Math.sqrt((state.x1 - state.x2) ** 2 + (state.y1 - state.y2) ** 2).toFixed(3)} m`
    });

    // 4. Charge Markers
    data.push({
      type: 'scatter',
      mode: 'markers+text',
      x: [state.x1, state.x2],
      y: [state.y1, state.y2],
      text: [`q₁ (${state.q1 >= 0 ? '+' : ''}${state.q1} nC)`, `q₂ (${state.q2 >= 0 ? '+' : ''}${state.q2} nC)`],
      textposition: ['top center', 'top center'],
      textfont: { family: 'Inter', size: 11, color: '#ffffff' },
      marker: {
        size: 18,
        color: [state.q1 >= 0 ? '#f43f5e' : '#06b6d4', state.q2 >= 0 ? '#f43f5e' : '#06b6d4'],
        line: { color: '#ffffff', width: 2.5 }
      },
      name: 'Cargas',
      showlegend: false
    });

    const layout = {
      ...plotlyCommonLayout,
      xaxis: {
        title: 'x (metros)',
        range: [-1.05, 1.05],
        zerolinecolor: '#475569',
        gridcolor: '#1e293b'
      },
      yaxis: {
        title: 'y (metros)',
        range: [-1.05, 1.05],
        zerolinecolor: '#475569',
        gridcolor: '#1e293b',
        scaleanchor: 'x',
        scaleratio: 1
      }
    };

    Plotly.react(elements.plotM2, data, layout, plotlyConfig);
  }

  // --- Module 3 Plot: Gauss's Law & Surface Invariance ---
  function renderModule3Plot(xArr, yArr, V) {
    if (!elements.plotM3) return;

    const data = [];

    // 1. Quiver vectors of electric field in background
    const quiver = buildQuiverVectors(17, [-0.95, 0.95], [-0.95, 0.95]);
    data.push({
      type: 'scatter',
      mode: 'lines',
      x: quiver.lineX,
      y: quiver.lineY,
      line: {
        color: 'rgba(56, 189, 248, 0.45)',
        width: 1.4
      },
      name: 'Líneas de Campo E',
      hoverinfo: 'none',
      showlegend: false
    });

    // 2. Gaussian Circle Perimeter & Enclosed Area
    const numCirclePts = 100;
    const circleX = [];
    const circleY = [];
    for (let k = 0; k <= numCirclePts; k++) {
      const theta = (k / numCirclePts) * 2 * Math.PI;
      circleX.push(state.radius * Math.cos(theta));
      circleY.push(state.radius * Math.sin(theta));
    }

    data.push({
      type: 'scatter',
      mode: 'lines',
      x: circleX,
      y: circleY,
      fill: 'toself',
      fillcolor: 'rgba(16, 185, 129, 0.12)',
      line: {
        color: '#10b981',
        width: 3,
        dash: 'dash'
      },
      name: `Superficie Gaussiana (R = ${state.radius.toFixed(2)} m)`,
      hoverinfo: 'name'
    });

    // 3. Normal vectors (dA) outward arrows along Gaussian boundary
    const normalCount = 12;
    const normLineX = [];
    const normLineY = [];
    const normLen = 0.08;
    for (let k = 0; k < normalCount; k++) {
      const theta = (k / normalCount) * 2 * Math.PI;
      const x0 = state.radius * Math.cos(theta);
      const y0 = state.radius * Math.sin(theta);
      const x1 = (state.radius + normLen) * Math.cos(theta);
      const y1 = (state.radius + normLen) * Math.sin(theta);

      normLineX.push(x0, x1, null);
      normLineY.push(y0, y1, null);
    }

    data.push({
      type: 'scatter',
      mode: 'lines',
      x: normLineX,
      y: normLineY,
      line: {
        color: '#a7f3d0',
        width: 2
      },
      name: 'Vectores Normales dA',
      hoverinfo: 'none'
    });

    // 4. Charge Markers with Enclosed status
    const d1 = Math.sqrt(state.x1 ** 2 + state.y1 ** 2);
    const d2 = Math.sqrt(state.x2 ** 2 + state.y2 ** 2);
    const q1Inside = d1 <= state.radius;
    const q2Inside = d2 <= state.radius;

    data.push({
      type: 'scatter',
      mode: 'markers+text',
      x: [state.x1, state.x2],
      y: [state.y1, state.y2],
      text: [
        `q₁ (${state.q1 >= 0 ? '+' : ''}${state.q1} nC)\n[${q1Inside ? 'DENTRO' : 'FUERA'}]`,
        `q₂ (${state.q2 >= 0 ? '+' : ''}${state.q2} nC)\n[${q2Inside ? 'DENTRO' : 'FUERA'}]`
      ],
      textposition: 'top center',
      textfont: { family: 'Inter', size: 10, color: '#f8fafc' },
      marker: {
        size: 18,
        color: [
          q1Inside ? '#f43f5e' : 'rgba(244,63,94,0.4)',
          q2Inside ? '#06b6d4' : 'rgba(6,182,212,0.4)'
        ],
        line: {
          color: [q1Inside ? '#ffffff' : '#f43f5e', q2Inside ? '#ffffff' : '#06b6d4'],
          width: 2.5
        }
      },
      name: 'Cargas Eléctricas',
      showlegend: false
    });

    const layout = {
      ...plotlyCommonLayout,
      xaxis: {
        title: 'x (metros)',
        range: [-1.05, 1.05],
        zerolinecolor: '#475569',
        gridcolor: '#1e293b'
      },
      yaxis: {
        title: 'y (metros)',
        range: [-1.05, 1.05],
        zerolinecolor: '#475569',
        gridcolor: '#1e293b',
        scaleanchor: 'x',
        scaleratio: 1
      }
    };

    Plotly.react(elements.plotM3, data, layout, plotlyConfig);
  }

  // --- CSV Export & Import Features ---

  /**
   * Generates CSV string of the grid mesh with metadata
   */
  function generateCSVData() {
    const N = state.resolution;
    const xArr = linspace(-1.0, 1.0, N);
    const yArr = linspace(-1.0, 1.0, N);
    const { V } = computePotentialGrid(xArr, yArr, state.q1, state.x1, state.y1, state.q2, state.x2, state.y2);

    let csv = `# GradVolt Web 3D - Sears & Zemansky Electrostatic Mesh Export\n`;
    csv += `# Parameters: q1=${state.q1},x1=${state.x1},y1=${state.y1},q2=${state.q2},x2=${state.x2},y2=${state.y2},R=${state.radius},N=${N}\n`;
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
    a.download = `gradvolt_potencial_malla_${state.resolution}x${state.resolution}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);

    showToast('Malla exportada con éxito en formato CSV', 'success');
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

        // Check for metadata line
        let metaParsed = false;
        for (let i = 0; i < Math.min(5, lines.length); i++) {
          if (lines[i].startsWith('#') && lines[i].includes('q1=')) {
            const parts = lines[i].replace('#', '').replace('Parameters:', '').trim().split(',');
            parts.forEach(p => {
              const [k, v] = p.split('=').map(s => s.trim());
              const val = parseFloat(v);
              if (!isNaN(val)) {
                if (k === 'q1') state.q1 = val;
                if (k === 'x1') state.x1 = val;
                if (k === 'y1') state.y1 = val;
                if (k === 'q2') state.q2 = val;
                if (k === 'x2') state.x2 = val;
                if (k === 'y2') state.y2 = val;
                if (k === 'R') state.radius = val;
                if (k === 'N') state.resolution = Math.round(val);
              }
            });
            metaParsed = true;
            break;
          }
        }

        // Parse data rows
        const xSet = new Set();
        const ySet = new Set();
        const rows = [];

        lines.forEach(line => {
          if (line.startsWith('#') || line.toLowerCase().startsWith('x,')) return;
          const cols = line.split(',').map(s => parseFloat(s.trim()));
          if (cols.length >= 3 && !isNaN(cols[0]) && !isNaN(cols[1]) && !isNaN(cols[2])) {
            xSet.add(cols[0]);
            ySet.add(cols[1]);
            rows.push({ x: cols[0], y: cols[1], v: cols[2] });
          }
        });

        if (rows.length === 0) {
          throw new Error('No se encontraron columnas válidas (x, y, V) en el CSV.');
        }

        state.importedFromCsv = false; // Always re-synced with mathematical model or metadata
        elements.selectResolution.value = String(state.resolution);

        scheduleUpdate();
        showToast(
          metaParsed 
            ? 'Simulación y parámetros restaurados desde CSV con éxito.' 
            : `Datos de potencial importados (${rows.length} puntos procesados).`, 
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

    elements.tablePointsCount.textContent = `(${xArr.length} × ${yArr.length} = ${xArr.length * yArr.length} celdas)`;

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
        // Reto 1: Ley Inverso Cuadrado (Carga única en el origen)
        state.q1 = 8.0;
        state.x1 = 0.00;
        state.y1 = 0.00;
        state.q2 = 0.0;
        state.x2 = 0.50;
        state.y2 = 0.00;
        state.radius = 0.40;
        setActiveTab('module-1');
        showToast('Reto 1 Aplicado: Carga puntual única. Analice la caída radial 1/r (V) y 1/r² (E).', 'info');
        break;

      case 2:
        // Reto 2: Paradoja Esfera de Gauss
        state.q1 = 5.0;
        state.x1 = 0.00;
        state.y1 = 0.00; // Inside Gaussian sphere
        state.q2 = -5.0;
        state.x2 = 0.70;
        state.y2 = 0.70; // Outside: distance = sqrt(0.7^2 + 0.7^2) = ~0.99 m
        state.radius = 0.50;
        setActiveTab('module-3');
        showToast('Reto 2 Aplicado: ¡El flujo Φ_E no cambia al variar R mientras q₂ siga fuera!', 'info');
        break;

      case 3:
        // Reto 3: Dipolo Eléctrico (V=0, E ≠ 0 en el plano medio)
        state.q1 = 6.0;
        state.x1 = -0.40;
        state.y1 = 0.00;
        state.q2 = -6.0;
        state.x2 = 0.40;
        state.y2 = 0.00;
        state.radius = 0.60;
        setActiveTab('module-2');
        showToast('Reto 3 Aplicado: En el eje vertical x=0, V=0 V pero el campo E es horizontal y máximo.', 'info');
        break;
    }

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

    // Resize Plotly plots to match current container dimensions
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

    const icon = type === 'success' ? 'fa-circle-check text-emerald-400' : type === 'error' ? 'fa-triangle-exclamation text-rose-400' : 'fa-circle-info text-sky-400';
    toast.innerHTML = `<i class="fa-solid ${icon}"></i><span>${message}</span>`;

    elements.toastContainer.appendChild(toast);

    setTimeout(() => {
      toast.style.opacity = '0';
      toast.style.transform = 'translateY(10px)';
      setTimeout(() => toast.remove(), 300);
    }, 3800);
  }

  // --- Event Listeners Setup ---
  function setupEventListeners() {
    // Charge 1 Sliders
    elements.sliderQ1.addEventListener('input', e => {
      state.q1 = parseFloat(e.target.value);
      state.importedFromCsv = false;
      scheduleUpdate();
    });
    elements.sliderX1.addEventListener('input', e => {
      state.x1 = parseFloat(e.target.value);
      state.importedFromCsv = false;
      scheduleUpdate();
    });
    elements.sliderY1.addEventListener('input', e => {
      state.y1 = parseFloat(e.target.value);
      state.importedFromCsv = false;
      scheduleUpdate();
    });

    // Charge 2 Sliders
    elements.sliderQ2.addEventListener('input', e => {
      state.q2 = parseFloat(e.target.value);
      state.importedFromCsv = false;
      scheduleUpdate();
    });
    elements.sliderX2.addEventListener('input', e => {
      state.x2 = parseFloat(e.target.value);
      state.importedFromCsv = false;
      scheduleUpdate();
    });
    elements.sliderY2.addEventListener('input', e => {
      state.y2 = parseFloat(e.target.value);
      state.importedFromCsv = false;
      scheduleUpdate();
    });

    // Gauss Radius Slider
    elements.sliderRadius.addEventListener('input', e => {
      state.radius = parseFloat(e.target.value);
      scheduleUpdate();
    });

    // Resolution selector
    elements.selectResolution.addEventListener('change', e => {
      state.resolution = parseInt(e.target.value, 10);
      scheduleUpdate();
    });

    // Reset defaults
    elements.btnResetDefaults.addEventListener('click', () => {
      state.q1 = 5.0;
      state.x1 = -0.40;
      state.y1 = 0.00;
      state.q2 = -5.0;
      state.x2 = 0.40;
      state.y2 = 0.00;
      state.radius = 0.50;
      state.resolution = 41;
      state.importedFromCsv = false;
      elements.selectResolution.value = '41';
      scheduleUpdate();
      showToast('Parámetros restablecidos a los valores por defecto.', 'info');
    });

    // Module 1 View Mode toggles
    elements.btnView3D.addEventListener('click', () => {
      state.m1ViewMode = '3d';
      elements.btnView3D.className = 'px-3 py-1 rounded-lg text-xs font-medium bg-brand-600 text-white shadow transition-all flex items-center gap-1.5';
      elements.btnViewContour.className = 'px-3 py-1 rounded-lg text-xs font-medium text-slate-400 hover:text-white transition-all flex items-center gap-1.5';
      scheduleUpdate();
    });

    elements.btnViewContour.addEventListener('click', () => {
      state.m1ViewMode = 'contour';
      elements.btnViewContour.className = 'px-3 py-1 rounded-lg text-xs font-medium bg-brand-600 text-white shadow transition-all flex items-center gap-1.5';
      elements.btnView3D.className = 'px-3 py-1 rounded-lg text-xs font-medium text-slate-400 hover:text-white transition-all flex items-center gap-1.5';
      scheduleUpdate();
    });

    elements.selectColorscale.addEventListener('change', e => {
      state.colorscale = e.target.value;
      scheduleUpdate();
    });

    // Module 1 Table Toggle
    elements.toggleTableBtn.addEventListener('click', () => {
      const isHidden = elements.tableContainer.classList.contains('hidden');
      elements.tableContainer.classList.toggle('hidden');
      elements.tableChevron.style.transform = isHidden ? 'rotate(180deg)' : 'rotate(0deg)';
    });

    // Module 2 Vector Options
    elements.checkShowEquipotential.addEventListener('change', e => {
      state.showEquipotential = e.target.checked;
      scheduleUpdate();
    });

    elements.checkShowVectors.addEventListener('change', e => {
      state.showVectors = e.target.checked;
      scheduleUpdate();
    });

    elements.selectVectorDensity.addEventListener('change', e => {
      state.vectorDensity = parseInt(e.target.value, 10);
      scheduleUpdate();
    });

    // Tab buttons
    document.querySelectorAll('.tab-btn').forEach(btn => {
      btn.addEventListener('click', () => {
        const tab = btn.getAttribute('data-tab');
        setActiveTab(tab);
      });
    });

    // Quick Scenarios
    elements.btnScenario1.addEventListener('click', () => applyScenario(1));
    elements.btnScenario2.addEventListener('click', () => applyScenario(2));
    elements.btnScenario3.addEventListener('click', () => applyScenario(3));

    // CSV Export & Import
    elements.btnExportCsv.addEventListener('click', exportCSV);
    if (elements.btnExportCsvSub) elements.btnExportCsvSub.addEventListener('click', exportCSV);

    elements.inputImportCsv.addEventListener('change', e => {
      const file = e.target.files[0];
      if (file) {
        importCSV(file);
        e.target.value = ''; // Reset input to allow re-uploading same file
      }
    });

    // Theory Modal
    elements.btnOpenTheory.addEventListener('click', () => {
      elements.modalTheory.classList.remove('hidden');
      renderMathInElement();
    });
    const closeModal = () => elements.modalTheory.classList.add('hidden');
    elements.btnCloseTheory.addEventListener('click', closeModal);
    elements.btnCloseTheory2.addEventListener('click', closeModal);
    elements.modalTheory.addEventListener('click', e => {
      if (e.target === elements.modalTheory) closeModal();
    });

    // Window Resize
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
    setupEventListeners();
    scheduleUpdate();

    // Trigger KaTeX after DOM is fully parsed
    setTimeout(() => {
      renderMathInElement();
    }, 400);
  }

  // Start on DOM ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }

})();
