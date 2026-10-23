/**
 * ============================================================================
 * GradVolt Web 3D - Suite de Pruebas Unitarias y de Integración DOM
 * Basado en Física Universitaria con Física Moderna de Sears & Zemansky
 * Capítulos 21 (Coulomb y Campo), 22 (Ley de Gauss) y 23 (Potencial)
 * 
 * Implementado con node:test y node:assert/strict (sin dependencias externas)
 * ============================================================================
 */

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');

// --- Constantes Físicas Fundamentales (SI) ---
const K_COULOMB = 8.9875517923e9; // Constante electrostática k = 1 / (4πε₀) [N·m²/C²]
const EPSILON_0 = 8.8541878176e-12; // Permitividad del vacío ε₀ [C²/(N·m²)]
const EPS_SMOOTH = 0.05; // Parámetro de suavizado (m) de la simulación

/**
 * Evalúa el potencial electrostático escalar V(x, y, z)
 * V(x, y, z) = sum_i [ k * q_i / sqrt((x-x_i)^2 + (y-y_i)^2 + z^2 + eps^2) ]
 */
function evaluatePotentialAt(charges, x, y, z = 0, eps = 0) {
  let vTotal = 0;
  for (const c of charges) {
    const qSI = c.q * 1e-9;
    const zi = c.z || 0;
    const r = Math.sqrt((x - c.x) ** 2 + (y - c.y) ** 2 + (z - zi) ** 2 + eps * eps);
    vTotal += K_COULOMB * (qSI / r);
  }
  return vTotal;
}

/**
 * Evalúa el vector campo eléctrico E(x, y, z) por superposición
 * E = sum_i [ k * q_i * (r - r_i) / |r - r_i|^3 ]
 */
function evaluateElectricFieldAt(charges, x, y, z = 0, eps = 0) {
  let ex = 0, ey = 0, ez = 0;
  for (const c of charges) {
    const qSI = c.q * 1e-9;
    const zi = c.z || 0;
    const dx = x - c.x;
    const dy = y - c.y;
    const dz = z - zi;
    const r2 = dx * dx + dy * dy + dz * dz + eps * eps;
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
 * Evalúa la Ley de Gauss para una superficie esférica centrada en el origen
 * Q_enc = sum_{r_i <= R} q_i
 * Phi_E = Q_enc / ε₀
 */
function computeGauss(charges, R) {
  let qEnc = 0;
  const inside = [];
  const outside = [];

  for (const c of charges) {
    const dist = Math.sqrt(c.x * c.x + c.y * c.y + (c.z || 0) ** 2);
    if (dist <= R) {
      qEnc += c.q;
      inside.push(c);
    } else {
      outside.push(c);
    }
  }

  const qEncSI = qEnc * 1e-9;
  const flux = qEncSI / EPSILON_0;
  const area = 4 * Math.PI * (R * R);
  const eAvg = Math.abs(flux) / area;

  return { qEnc, qEncSI, flux, area, eAvg, inside, outside };
}

// ============================================================================
// BLOQUE DE PRUEBAS
// ============================================================================

test('Prueba 1: Reto 1 - Ley del Inverso del Cuadrado (Sears & Zemansky Caps. 21 y 23)', () => {
  // Configuración de 1 carga puntual única en el origen (q = +8.0 nC)
  const singleCharge = [
    { id: 1, q: 8.0, x: 0.00, y: 0.00, z: 0.00 }
  ];

  const r1 = 0.25; // metros
  const r2 = 0.50; // metros (duplicamos la distancia: r2 = 2 * r1)

  // 1. Evaluación del Potencial Eléctrico Escalar V: decae como 1/r
  const V1 = evaluatePotentialAt(singleCharge, r1, 0, 0, 0);
  const V2 = evaluatePotentialAt(singleCharge, r2, 0, 0, 0);

  // Al duplicar r, V debe reducirse a la mitad exacta (V2 = V1 / 2)
  const ratioV = V2 / V1;
  assert.ok(
    Math.abs(ratioV - 0.5) < 1e-7,
    `El potencial debe decaer como 1/r: V(0.50) / V(0.25) = ${ratioV}, esperado 0.5`
  );
  assert.ok(
    Math.abs(V2 - V1 / 2) < 1e-4,
    `V2 (${V2} V) debe ser exactamente la mitad de V1 (${V1} V)`
  );

  // 2. Evaluación del Campo Eléctrico Vectorial |E|: decae como 1/r²
  const E1 = evaluateElectricFieldAt(singleCharge, r1, 0, 0, 0);
  const E2 = evaluateElectricFieldAt(singleCharge, r2, 0, 0, 0);

  // Al duplicar r, |E| debe reducirse a la cuarta parte exacta (|E2| = |E1| / 4)
  const ratioE = E2.mag / E1.mag;
  assert.ok(
    Math.abs(ratioE - 0.25) < 1e-7,
    `El campo eléctrico debe decaer como 1/r²: |E(0.50)| / |E(0.25)| = ${ratioE}, esperado 0.25`
  );
  assert.ok(
    Math.abs(E2.mag - E1.mag / 4) < 1e-3,
    `|E2| (${E2.mag} N/C) debe ser exactamente un cuarto de |E1| (${E1.mag} N/C)`
  );

  // Verificamos que la dirección del campo sea estrictamente radial hacia afuera (+x)
  assert.ok(E1.ex > 0 && Math.abs(E1.ey) < 1e-10, 'El campo debe ser netamente radial en +x');
});

test('Prueba 2: Reto 2 - Invarianza de la Ley de Gauss (Sears & Zemansky Cap. 22)', () => {
  // Configuración del Reto 2:
  // q1 = +5.0 nC en el origen (0, 0)
  // q2 = -5.0 nC en (0.70, 0.70) -> distancia d2 = sqrt(0.70^2 + 0.70^2) ≈ 0.9899 m
  const charges = [
    { id: 1, q: 5.0, x: 0.00, y: 0.00, z: 0.00 },
    { id: 2, q: -5.0, x: 0.70, y: 0.70, z: 0.00 }
  ];

  // Caso A: Radio R1 = 0.30 m
  const g1 = computeGauss(charges, 0.30);
  assert.equal(g1.inside.length, 1, 'En R = 0.30 m solo debe estar dentro q1');
  assert.equal(g1.outside.length, 1, 'q2 debe estar fuera de R = 0.30 m');
  assert.equal(g1.qEnc, 5.0, 'Q_enc en R = 0.30 m debe ser +5.0 nC');

  // Caso B: Radio R2 = 0.60 m (mayor radio, pero q2 sigue afuera ya que d2 ≈ 0.99 m > 0.60 m)
  const g2 = computeGauss(charges, 0.60);
  assert.equal(g2.inside.length, 1, 'En R = 0.60 m solo debe estar dentro q1');
  assert.equal(g2.qEnc, 5.0, 'Q_enc en R = 0.60 m debe seguir siendo +5.0 nC');

  // Demostración de invarianza de flujo Phi_E = Q_enc / ε₀
  assert.equal(
    g1.flux,
    g2.flux,
    `El flujo debe ser estrictamente idéntico (Invarianza de Gauss): Phi1 = ${g1.flux}, Phi2 = ${g2.flux}`
  );
  assert.ok(g1.flux > 0, 'El flujo neto debe ser positivo ya que Q_enc > 0');

  // Caso C: Radio R3 = 1.10 m (encierra ambas cargas, ya que d2 ≈ 0.99 m < 1.10 m)
  const g3 = computeGauss(charges, 1.10);
  assert.equal(g3.inside.length, 2, 'En R = 1.10 m deben estar encerradas ambas cargas (q1 y q2)');
  assert.equal(g3.qEnc, 0.0, 'La carga neta encerrada Q_enc debe ser exactamente 0 nC (5.0 - 5.0)');
  assert.equal(g3.flux, 0.0, 'El flujo eléctrico total debe ser rigurosamente 0 N·m²/C');
});

test('Prueba 3: Reto 3 - Dipolo Eléctrico Simétrico (Sears & Zemansky Caps. 21 y 23)', () => {
  // Configuración del Reto 3: Dipolo simétrico sobre el eje x
  // q1 = +6.0 nC en (-0.40, 0)
  // q2 = -6.0 nC en (+0.40, 0)
  const dipole = [
    { id: 1, q: 6.0, x: -0.40, y: 0.00, z: 0.00 },
    { id: 2, q: -6.0, x: 0.40, y: 0.00, z: 0.00 }
  ];

  // 1. En el origen (0, 0), equidistante a ambas cargas:
  // V(0,0) = k*(+q)/d + k*(-q)/d = 0 V
  const V_origin = evaluatePotentialAt(dipole, 0.0, 0.0, 0, 0);
  assert.equal(
    V_origin,
    0.0,
    `El potencial en el plano medio simétrico debe ser exactamente 0 V, obtenido: ${V_origin}`
  );

  // 2. El campo eléctrico en el origen NO es nulo; las contribuciones se suman constructivamente en +x:
  // E_total = E1 (+x) + E2 (+x) = 2 * k * q / d^2 en dirección +x
  const E_origin = evaluateElectricFieldAt(dipole, 0.0, 0.0, 0, 0);
  
  assert.ok(
    E_origin.ex > 0,
    `El campo eléctrico horizontal debe apuntar hacia la carga negativa (+x): Ex = ${E_origin.ex}`
  );
  assert.equal(
    Math.abs(E_origin.ey),
    0.0,
    'Por simetría del dipolo, la componente vertical Ey debe ser estrictamente cero'
  );
  assert.ok(
    E_origin.mag > 0,
    `La magnitud del campo eléctrico debe ser distinta de cero (|E| = ${E_origin.mag} N/C)`
  );

  // Verificación cuantitativa teórica del dipolo: E_teorico = 2 * k * q / d^2
  const d = 0.40;
  const E_expected = 2 * (K_COULOMB * (6.0 * 1e-9)) / (d * d);
  assert.ok(
    Math.abs(E_origin.ex - E_expected) < 1e-4,
    `El campo medido (${E_origin.ex}) debe coincidir con el valor analítico (${E_expected} N/C)`
  );
});

test('Prueba 4: Superposición Vectorial y Escalar con 10 Cargas Simultáneas', () => {
  // Arreglo completo con las 10 cargas puntuales activas
  const tenCharges = [
    { id: 1, q: 5.0, x: -0.40, y: 0.00 },
    { id: 2, q: -5.0, x: 0.40, y: 0.00 },
    { id: 3, q: 3.5, x: 0.00, y: 0.35 },
    { id: 4, q: -2.0, x: 0.25, y: -0.45 },
    { id: 5, q: 4.0, x: -0.30, y: 0.60 },
    { id: 6, q: -1.5, x: 0.50, y: 0.20 },
    { id: 7, q: 6.0, x: -0.60, y: -0.30 },
    { id: 8, q: -4.5, x: 0.15, y: 0.70 },
    { id: 9, q: 2.0, x: -0.20, y: -0.70 },
    { id: 10, q: -3.0, x: 0.80, y: -0.25 }
  ];

  assert.equal(tenCharges.length, 10, 'Deben haber exactamente 10 cargas simultáneas');

  // Evaluar en puntos de prueba del espacio
  const testPoints = [
    { x: 0.00, y: 0.00 },
    { x: 0.15, y: -0.20 },
    { x: -0.50, y: 0.50 },
    { x: 0.75, y: 0.75 }
  ];

  for (const pt of testPoints) {
    const V = evaluatePotentialAt(tenCharges, pt.x, pt.y, 0, EPS_SMOOTH);
    const E = evaluateElectricFieldAt(tenCharges, pt.x, pt.y, 0, EPS_SMOOTH);

    assert.ok(Number.isFinite(V), `El potencial en (${pt.x}, ${pt.y}) debe ser finito`);
    assert.ok(Number.isFinite(E.ex), `Ex en (${pt.x}, ${pt.y}) debe ser finito`);
    assert.ok(Number.isFinite(E.ey), `Ey en (${pt.x}, ${pt.y}) debe ser finito`);
    assert.ok(Number.isFinite(E.mag), `|E| en (${pt.x}, ${pt.y}) debe ser finito`);
    assert.ok(E.mag >= 0, `|E| en (${pt.x}, ${pt.y}) debe ser no negativo`);
  }

  // Evaluar Ley de Gauss con 10 cargas para radio R = 0.50 m
  const gauss = computeGauss(tenCharges, 0.50);
  assert.ok(Number.isFinite(gauss.qEnc), 'Q_enc con 10 cargas debe ser un número finito');
  assert.ok(Number.isFinite(gauss.flux), 'El flujo Phi_E con 10 cargas debe ser un número finito');
  assert.ok(Number.isFinite(gauss.eAvg), 'El campo medio superficial debe ser un número finito');
  assert.equal(
    gauss.inside.length + gauss.outside.length,
    10,
    'La suma de cargas dentro y fuera debe ser exactamente 10'
  );
});

test('Prueba 5: Integridad Estructural de Interfaz HTML y Lógica JS (Plotly Purge & DOM)', () => {
  const htmlPath = path.join(__dirname, '..', 'index.html');
  const jsPath = path.join(__dirname, '..', 'app.js');

  assert.ok(fs.existsSync(htmlPath), 'index.html debe existir en la raíz del proyecto');
  assert.ok(fs.existsSync(jsPath), 'app.js debe existir en la raíz del proyecto');

  const htmlContent = fs.readFileSync(htmlPath, 'utf8');
  const jsContent = fs.readFileSync(jsPath, 'utf8');

  // 1. Verificación de IDs esenciales en index.html
  const requiredIds = [
    'plot-module-1',
    'plot-module-2',
    'plot-module-3',
    'btn-scenario-1',
    'btn-scenario-2',
    'btn-scenario-3',
    'btn-export-csv',
    'input-import-csv',
    'btn-add-charge',
    'charges-list-container'
  ];

  for (const id of requiredIds) {
    const idPattern = new RegExp(`id=["']${id}["']`);
    assert.ok(
      idPattern.test(htmlContent),
      `index.html debe contener el elemento con id="${id}"`
    );
  }

  // 2. Verificación de purga de Plotly al cambiar entre 3D y 2D en app.js
  assert.ok(
    jsContent.includes('Plotly.purge'),
    'app.js debe invocar Plotly.purge para evitar colisiones entre cámara 3D y ejes 2D'
  );

  // 3. Verificación de que el error crítico de plotly_click y variables sobrantes NO están presentes
  assert.ok(
    !jsContent.includes("plotEl.on('plotly_click'"),
    'app.js NO debe registrar eventos plotEl.on sobre elementos de Plotly sin inicializar'
  );
  assert.ok(
    !jsContent.includes('clickToPlaceMode'),
    'app.js no debe contener clickToPlaceMode tras la limpieza de la interfaz'
  );
  assert.ok(
    !htmlContent.includes('btn-mode-2d') && !htmlContent.includes('btn-mode-3d'),
    'index.html no debe contener los botones eliminados de dimensión 3D'
  );
});
