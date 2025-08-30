// ---------- Utility helpers ----------
const R_UNIVERSAL = 8.31446261815324e3; // J / (kmol K)  (also equals 8.31446261815324 J/(mol K) * 1000)

/** format numbers nicely or show '—' */
const fmt = (x, digits = 3) =>
  Number.isFinite(x) ? Number(x).toLocaleString(undefined, { maximumFractionDigits: digits }) : '—';

/** safe parse float from selector; returns NaN if not present or empty */
function getFloat(sel) {
  const el = document.querySelector(sel);
  if (!el) return NaN;
  const v = parseFloat(el.value);
  return Number.isFinite(v) ? v : NaN;
}

/** convert MPa->Pa if value looks like MPa (we assume user enters MPa as in UI) */
function MPaToPa(x) {
  if (!Number.isFinite(x)) return NaN;
  return x * 1e6;
}

// ---------- Isentropic area relation ----------
/*
  A/A* = (1/M) * [ (2/(gamma+1)) * (1 + (gamma-1)/2 * M^2) ]^{(gamma+1)/(2(gamma-1))}
*/
function areaRatioFromMach(M, gamma) {
  const g = gamma;
  const term = (2 / (g + 1)) * (1 + (g - 1) / 2 * M * M);
  const pow = (g + 1) / (2 * (g - 1));
  return (1 / M) * Math.pow(term, pow);
}

/** solve Mach for a given area ratio:
    - prefer the supersonic solution for areaRatio >= 1
    - numeric bisection (robust)
*/
function machFromAreaRatio(targetAR, gamma) {
  if (!Number.isFinite(targetAR) || !Number.isFinite(gamma)) return NaN;
  // Handle trivial:
  if (targetAR <= 1) {
    // subsonic branch (M < 1): search (1e-6, 0.9999)
    let lo = 1e-8, hi = 0.999999;
    for (let i = 0; i < 80; ++i) {
      const mid = 0.5 * (lo + hi);
      const f = areaRatioFromMach(mid, gamma) - targetAR;
      if (Math.abs(f) < 1e-9) return mid;
      // monotonic: area increases as M decreases below 1? for subsonic area decreases with M? safe bisection:
      if (areaRatioFromMach(lo, gamma) - targetAR > 0) { // adjust sign logic robustly
        if (f > 0) lo = mid; else hi = mid;
      } else {
        if (f > 0) hi = mid; else lo = mid;
      }
    }
    return 0.5 * (lo + hi);
  } else {
    // supersonic branch (M > 1): search (1.000001, large)
    let lo = 1.0000001, hi = 200; // hi large enough for rockets
    // Ensure function crosses sign:
    const f_lo = areaRatioFromMach(lo, gamma) - targetAR;
    const f_hi = areaRatioFromMach(hi, gamma) - targetAR;
    // if targetAR is very large we may need to increase hi; do it adaptively
    let hiAttempt = hi;
    let fhi = f_hi;
    let attempts = 0;
    while (f_lo * fhi > 0 && attempts < 10) {
      hiAttempt *= 2;
      fhi = areaRatioFromMach(hiAttempt, gamma) - targetAR;
      attempts++;
    }
    if (f_lo * fhi > 0) {
      // fallback: try simple Newton-like initial guess (not ideal). Return NaN if fail.
      return NaN;
    }
    hi = hiAttempt;
    // bisection:
    for (let i = 0; i < 80; ++i) {
      const mid = 0.5 * (lo + hi);
      const fmid = areaRatioFromMach(mid, gamma) - targetAR;
      if (Math.abs(fmid) < 1e-10) return mid;
      if ((areaRatioFromMach(lo, gamma) - targetAR) * fmid <= 0) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    return 0.5 * (lo + hi);
  }
}

// ---------- Main performance calculation ----------
function calculateIsPerf() {
  // read mandatory inputs
  const expansionRatio = getFloat('#expansionRatio'); // Ae/At (-)
  const gamma = getFloat('#specificHeatRatio'); // gamma
  let Pc_MPa = getFloat('#chamberPressure'); // MPa
  const T0 = getFloat('#chamberTemperature'); // K

  // optional inputs
  let ambientPressure_MPa = getFloat('#ambientPressure'); // MPa (optional)
  const molarMass = getFloat('#molarMass'); // optional: in g/mol (or kg/kmol) — numeric value 28.97 for air
  const At = getFloat('#throatArea'); // optional throat area in m^2 (needed for mdot, thrust)

  // convert pressures to Pa (we assume user enters MPa, consistent with your UI)
  const P0 = Number.isFinite(Pc_MPa) ? MPaToPa(Pc_MPa) : NaN; // chamber stagnation pressure (Pa)
  const Pa = Number.isFinite(ambientPressure_MPa) ? MPaToPa(ambientPressure_MPa) : NaN;

  // Prepare results container (will create if missing)
  let results = document.querySelector('#perf-results');
  if (!results) {
    results = document.createElement('div');
    results.id = 'perf-results';
    results.style.marginTop = '1rem';
    results.style.padding = '1rem';
    results.style.background = 'rgba(255,255,255,0.02)';
    results.style.borderRadius = '8px';
    // append after the form or to body if not found
    const form = document.querySelector('#calc-form') || document.body;
    form.insertAdjacentElement('afterend', results);
  }
  results.innerHTML = '<h3>Performance results</h3>';

  // If gamma or expansion ratio missing -> nothing to compute
  if (!Number.isFinite(expansionRatio) || !Number.isFinite(gamma)) {
    results.innerHTML += `<p style="color:orange">Provide at least expansion ratio and specific heat ratio (gamma) to compute area/Mach relations.</p>`;
    return;
  }

  // 1) compute exit Mach from expansion ratio (prefer supersonic branch)
  const Me = machFromAreaRatio(expansionRatio, gamma);
  results.innerHTML += `<div>Exit Mach (M_e): <strong>${fmt(Me, 4)}</strong></div>`;

  // 2) stagnation ratios (Pe/P0, Te/T0)
  const Te_div_T0 = 1 / (1 + (gamma - 1) / 2 * Me * Me);
  const Pe_div_P0 = Math.pow(1 + (gamma - 1) / 2 * Me * Me, -gamma / (gamma - 1));
  results.innerHTML += `<div>Stagnation temperature ratio (T_e/T_0): <strong>${fmt(Te_div_T0,5)}</strong></div>`;
  results.innerHTML += `<div>Stagnation pressure ratio (P_e/P_0): <strong>${fmt(Pe_div_P0,5)}</strong></div>`;

  // 3) exit static temperature & pressure if total chamber values given
  const Te = Number.isFinite(T0) ? T0 * Te_div_T0 : NaN;
  const Pe = Number.isFinite(P0) ? P0 * Pe_div_P0 : NaN;
  results.innerHTML += `<div>Exit temperature (T_e) [K]: <strong>${fmt(Te)}</strong></div>`;
  results.innerHTML += `<div>Exit pressure (P_e) [Pa]: <strong>${fmt(Pe)}</strong></div>`;

  // 4) exit velocity if R is known (need molarMass)
  let Rspec = NaN;
  if (Number.isFinite(molarMass)) {
    // molarMass: user input in "g/mol" or "kg/kmol" numeric; they are numerically equal (e.g., 28.97).
    // Rspec = R_universal / M (J / (kg K)), where M is in kg/kmol
    Rspec = R_UNIVERSAL / molarMass; // J/(kg K)
  }

  const Ve = Number.isFinite(Rspec) && Number.isFinite(Te) ? Me * Math.sqrt(gamma * Rspec * Te) : NaN;
  results.innerHTML += `<div>Exit velocity (V_e) [m/s]: <strong>${fmt(Ve)}</strong></div>`;

  // 5) Thrust coefficient Cf (isentropic ideal expression) - needs Pe, Pa, Pc and expansion ratio
  let Cf = NaN;
  if (Number.isFinite(P0)) {
    // need ambient pressure Pa to compute pressure term; if missing, assume vacuum Pa = 0 for term but show warning
    const Pa_for_Cf = Number.isFinite(Pa) ? Pa : 0;
    // term1: momentum term from expansion (isentropic)
    const term1 = Math.sqrt(
      (2 * Math.pow(gamma, 2)) / (gamma - 1) *
        Math.pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1)) *
        (1 - Math.pow(Pe_div_P0, (gamma - 1) / gamma))
    );
    const term2 = ((Pe - Pa_for_Cf) / P0) * expansionRatio; // pressure correction
    Cf = term1 + term2;
    results.innerHTML += `<div>Thrust coefficient (C_f): <strong>${fmt(Cf,4)}</strong> ${Number.isFinite(pa) ? '' : '<em>(pa assumed 0 for Cf)</em>'}</div>`;
  } else {
    results.innerHTML += `<div>Thrust coefficient (C_f): <strong>${fmt(Cf)}</strong> (needs chamber pressure P0)</div>`;
  }

  // 6) mass flow (choked throat) mdot if At and Rspec and T0 and P0 available
  let mdot = NaN;
  if (Number.isFinite(At) && Number.isFinite(Rspec) && Number.isFinite(T0) && Number.isFinite(P0)) {
    // choked mass flow formula for ideal gas at throat (M=1):
    // mdot = (P0 * At) * sqrt(gamma / (Rspec * T0)) * ( (2/(gamma+1))^((gamma+1)/(2*(gamma-1))) )
    const chokeFactor = Math.pow(2 / (gamma + 1), (gamma + 1) / (2 * (gamma - 1)));
    mdot = P0 * At * Math.sqrt(gamma / (Rspec * T0)) * chokeFactor;
    results.innerHTML += `<div>Mass flow (ṁ) [kg/s]: <strong>${fmt(mdot)}</strong></div>`;
  } else {
    results.innerHTML += `<div>Mass flow (ṁ): <strong>${fmt(mdot)}</strong> (needs throat area, molar mass, chamber P and T)</div>`;
  }

  // 7) Thrust and vacuum thrust if Cf and At and P0 known
  let thrust = NaN, thrustVac = NaN;
  if (Number.isFinite(Cf) && Number.isFinite(P0) && Number.isFinite(At)) {
    thrust = Cf * P0 * At;
    // vacuum thrust -> use Pa=0 in Cf formula: recompute Cf_vac (we already used Pa_for_Cf possibly zero)
    // For vacuum Cf we want explicit formula; simpler: recompute with Pa=0:
    const term1_v = Math.sqrt(
      (2 * Math.pow(gamma, 2)) / (gamma - 1) *
        Math.pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1)) *
        (1 - Math.pow(Pe_div_P0, (gamma - 1) / gamma))
    );
    const term2_v = ((Pe - 0) / P0) * expansionRatio;
    const Cf_vac = term1_v + term2_v;
    thrustVac = Cf_vac * P0 * At;
    results.innerHTML += `<div>Thrust (sea-level) [N]: <strong>${fmt(thrust)}</strong></div>`;
    results.innerHTML += `<div>Thrust (vacuum) [N]: <strong>${fmt(thrustVac)}</strong></div>`;
  } else {
    results.innerHTML += `<div>Thrust: <strong>${fmt(thrust)}</strong> (needs C_f, P0 and throat area)</div>`;
  }

  // 8) Isp if thrust and mdot available
  const g0 = 9.80665;
  let Isp = NaN;
  if (Number.isFinite(thrust) && Number.isFinite(mdot) && mdot > 0) {
    Isp = thrust / (mdot * g0);
    results.innerHTML += `<div>Specific impulse (I_sp) [s]: <strong>${fmt(Isp,3)}</strong></div>`;
  } else {
    results.innerHTML += `<div>Specific impulse (I_sp): <strong>${fmt(Isp)}</strong> (needs thrust & mass flow)</div>`;
  }

  // 9) Vacuum-specific outputs if possible
  // already provided thrustVac above; also compute vacuum Isp if possible:
  if (Number.isFinite(thrustVac) && Number.isFinite(mdot) && mdot > 0) {
    const Isp_vac = thrustVac / (mdot * g0);
    results.innerHTML += `<div>Vacuum I_sp [s]: <strong>${fmt(Isp_vac,3)}</strong></div>`;
  }

  // 10) final summary additions
  results.innerHTML += `<hr style="opacity:0.08">`;
  results.innerHTML += `<div>Inputs used:<pre style="white-space:pre-wrap">${JSON.stringify({
    expansionRatio, gamma, P0: Number.isFinite(P0) ? fmt(P0) + ' Pa' : 'n/a',
    T0: Number.isFinite(T0) ? fmt(T0)+' K' : 'n/a',
    molarMass: Number.isFinite(molarMass) ? molarMass+' (g/mol)' : 'n/a',
    throatArea: Number.isFinite(At) ? At+' m^2' : 'n/a',
    ambientPressure: Number.isFinite(Pa) ? fmt(Pa)+' Pa' : 'n/a'
  }, null, 2)}</pre></div>`;

  // Persist last values (only inputs)
  try {
    localStorage.setItem('calc-state', JSON.stringify({
      expansionRatio, specificHeatRatio: gamma,
      chamberPressure: Pc_MPa, chamberTemperature: T0,
      ambientPressure_MPa, molarMass, throatArea: At
    }));
  } catch {}
}
