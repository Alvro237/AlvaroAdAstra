// ---------- Utility helpers ----------
const R_UNIVERSAL = 8.31446261815324e3; // J / (kmol K)

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
  if (targetAR <= 1) {
    let lo = 1e-8, hi = 0.999999;
    for (let i = 0; i < 80; ++i) {
      const mid = 0.5 * (lo + hi);
      const f = areaRatioFromMach(mid, gamma) - targetAR;
      if (Math.abs(f) < 1e-9) return mid;
      if (areaRatioFromMach(lo, gamma) - targetAR > 0) {
        if (f > 0) lo = mid; else hi = mid;
      } else {
        if (f > 0) hi = mid; else lo = mid;
      }
    }
    return 0.5 * (lo + hi);
  } else {
    let lo = 1.0000001, hi = 200;
    const f_lo = areaRatioFromMach(lo, gamma) - targetAR;
    let hiAttempt = hi;
    let fhi = areaRatioFromMach(hiAttempt, gamma) - targetAR;
    let attempts = 0;
    while (f_lo * fhi > 0 && attempts < 10) {
      hiAttempt *= 2;
      fhi = areaRatioFromMach(hiAttempt, gamma) - targetAR;
      attempts++;
    }
    if (f_lo * fhi > 0) {
      return NaN;
    }
    hi = hiAttempt;
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

// --- Main performance calculation (FIXED) ---
function calculateIsPerf() {
  // read mandatory inputs
  const expansionRatio = getFloat('#expansionRatio'); // Ae/At (-)
  const gamma = getFloat('#specificHeatRatio'); // gamma
  let Pc_MPa = getFloat('#chamberPressure'); // MPa
  const T0 = getFloat('#chamberTemperature'); // K

  // optional inputs
  let ambientPressure_MPa = getFloat('#ambientPressure'); // MPa
  const molarMass = getFloat('#molarMass'); // optional: in g/mol
  const At = getFloat('#throatArea'); // optional throat area in m^2

  // convert pressures to Pa
  const P0 = Number.isFinite(Pc_MPa) ? MPaToPa(Pc_MPa) : NaN;
  const Pa = Number.isFinite(ambientPressure_MPa) ? MPaToPa(ambientPressure_MPa) : NaN;

  // Clear previous results in the table
  document.querySelectorAll('.results td').forEach(cell => cell.textContent = '—');

  // If gamma or expansion ratio missing -> nothing to compute
  if (!Number.isFinite(expansionRatio) || !Number.isFinite(gamma)) {
    alert('Please provide at least the expansion ratio and specific heat ratio.');
    return;
  }

  // 1) compute exit Mach from expansion ratio
  const Me = machFromAreaRatio(expansionRatio, gamma);

  // 2) stagnation ratios (Pe/P0, Te/T0)
  const Te_div_T0 = 1 / (1 + (gamma - 1) / 2 * Me * Me);
  const Pe_div_P0 = Math.pow(1 + (gamma - 1) / 2 * Me * Me, -gamma / (gamma - 1));

  // 3) exit static temperature & pressure if total chamber values given
  const Te = Number.isFinite(T0) ? T0 * Te_div_T0 : NaN;
  const Pe = Number.isFinite(P0) ? P0 * Pe_div_P0 : NaN;

  // 4) exit velocity if R is known
  let Rspec = NaN;
  if (Number.isFinite(molarMass)) {
    Rspec = R_UNIVERSAL / molarMass; // J/(kg K)
  }
  const Ve = Number.isFinite(Rspec) && Number.isFinite(Te) ? Me * Math.sqrt(gamma * Rspec * Te) : NaN;

  // 5) Thrust coefficient Cf
  let Cf = NaN;
  const Pe_div_P0_term = Math.pow(Pe_div_P0, (gamma - 1) / gamma);
  if (Number.isFinite(P0)) {
    const Pa_for_Cf = Number.isFinite(Pa) ? Pa : 0;
    const term1 = Math.sqrt(
      (2 * Math.pow(gamma, 2)) / (gamma - 1) *
      Math.pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1)) *
      (1 - Pe_div_P0_term)
    );
    const term2 = ((Pe - Pa_for_Cf) / P0) * expansionRatio;
    Cf = term1 + term2;
  }

  // 6) mass flow
  let mdot = NaN;
  if (Number.isFinite(At) && Number.isFinite(Rspec) && Number.isFinite(T0) && Number.isFinite(P0)) {
    const chokeFactor = Math.pow(2 / (gamma + 1), (gamma + 1) / (2 * (gamma - 1)));
    mdot = P0 * At * Math.sqrt(gamma / (Rspec * T0)) * chokeFactor;
  }

  // 7) Thrust and vacuum thrust
  let thrust = NaN, thrustVac = NaN;
  if (Number.isFinite(Cf) && Number.isFinite(P0) && Number.isFinite(At)) {
    thrust = Cf * P0 * At;
    const term1_v = Math.sqrt(
      (2 * Math.pow(gamma, 2)) / (gamma - 1) *
      Math.pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1)) *
      (1 - Pe_div_P0_term)
    );
    const term2_v = (Pe / P0) * expansionRatio;
    const Cf_vac = term1_v + term2_v;
    thrustVac = Cf_vac * P0 * At;
  }

  // 8) Isp
  const g0 = 9.80665;
  let Isp = NaN;
  if (Number.isFinite(thrust) && Number.isFinite(mdot) && mdot > 0) {
    Isp = thrust / (mdot * g0);
  }

  // 9) Populate the HTML table
  document.getElementById('vExit').textContent = fmt(Ve);
  document.getElementById('pExit').textContent = fmt(Pe);
  document.getElementById('tExit').textContent = fmt(Te);
  document.getElementById('machExit').textContent = fmt(Me, 4);
  document.getElementById('pRatio').textContent = fmt(Pe_div_P0, 5);
  document.getElementById('tRatio').textContent = fmt(Te_div_T0, 5);
  document.getElementById('cf').textContent = fmt(Cf, 4);
  document.getElementById('thrust').textContent = fmt(thrust);
  document.getElementById('thrustVac').textContent = fmt(thrustVac);
  document.getElementById('isp').textContent = fmt(Isp, 3);
  document.getElementById('mdot').textContent = fmt(mdot);

  // Optional: You can still show a summary if you want, but the main results are in the table now.
  // The 'perf-results' div from your original code is not in your HTML, so this part is commented out.
  /*
  const summaryDiv = document.querySelector('.results');
  if (summaryDiv) {
      summaryDiv.innerHTML += `<div>Inputs used: <pre>${JSON.stringify({
          expansionRatio, gamma, P0, T0, molarMass, throatArea: At, ambientPressure: Pa
      }, null, 2)}</pre></div>`;
  }
  */
}

document.getElementById("compute").addEventListener("click", calculateIsPerf);



  // ---------- Isentropic relations ----------
    function ratiosFromMach(M, gamma) {
      const g = gamma;
      const T0T = 1 + ((g - 1) / 2) * M * M;
      const p0p = Math.pow(T0T, g / (g - 1));
      const rho0rho = Math.pow(T0T, 1 / (g - 1));
      const AAstar = areaRatioFromMach(M, g);
      const TTstar = 1 / T0T * (1 + (g - 1) / 2);
      const ppstar = (1 / p0p) / Math.pow(2 / (g + 1), g / (g - 1));
      const rhorhostar = 1 / rho0rho / Math.pow(2 / (g + 1), 1 / (g - 1));
      return {T0T, p0p, rho0rho, AAstar, TTstar, ppstar, rhorhostar};
    }

    function areaRatioFromMach(M, g) {
      const term = (2 / (g + 1)) * (1 + (g - 1) / 2 * M * M);
      const pow = (g + 1) / (2 * (g - 1));
      return (1 / M) * Math.pow(term, pow);
    }

    function machFromT0T(val, g) {
      return Math.sqrt(2 * (val - 1) / (g - 1));
    }
    function machFromP0P(val, g) {
      const T0T = Math.pow(val, (g - 1) / g);
      return machFromT0T(T0T, g);
    }
    function machFromRho0Rho(val, g) {
      const T0T = Math.pow(val, g - 1);
      return machFromT0T(T0T, g);
    }
    function machFromAreaRatio(val, g) {
      return machFromAreaRatioHelper(val, g); // reusing your bisection solver if desired
    }

    // Minimal bisection solver for A/A*
    function machFromAreaRatioHelper(targetAR, g) {
      let lo = 1e-6, hi = 20;
      for (let i = 0; i < 100; i++) {
        const mid = 0.5 * (lo + hi);
        const f = areaRatioFromMach(mid, g) - targetAR;
        if (Math.abs(f) < 1e-6) return mid;
        if (f > 0) hi = mid; else lo = mid;
      }
      return 0.5 * (lo + hi);
    }

    // ---------- Main function ----------
    function fillIsentropicTable() {
      const g = parseFloat(document.getElementById('gamma').value);
      const param = document.getElementById('knownParam').value;
      const val = parseFloat(document.getElementById('paramValue').value);

      if (!Number.isFinite(g) || !Number.isFinite(val)) {
        alert("Please enter valid inputs.");
        return;
      }

      let M;
      if (param === "Ma") M = val;
      if (param === "T0T") M = machFromT0T(val, g);
      if (param === "p0p") M = machFromP0P(val, g);
      if (param === "rho0rho") M = machFromRho0Rho(val, g);
      if (param === "AAstar") M = machFromAreaRatioHelper(val, g);
      // TT*, p/p*, ρ/ρ* inversion could also be added, but usually we define Mach from T0/T, p0/p, etc.

      if (!Number.isFinite(M)) {
        alert("Could not determine Mach number from given input.");
        return;
      }

      const {T0T, p0p, rho0rho, AAstar, TTstar, ppstar, rhorhostar} = ratiosFromMach(M, g);

      document.getElementById('res-gamma').textContent = fmt(g);
      document.getElementById('res-Ma').textContent = fmt(M);
      document.getElementById('res-T0T').textContent = fmt(T0T);
      document.getElementById('res-p0p').textContent = fmt(p0p);
      document.getElementById('res-rho0rho').textContent = fmt(rho0rho);
      document.getElementById('res-AAstar').textContent = fmt(AAstar);
      document.getElementById('res-TTstar').textContent = fmt(TTstar);
      document.getElementById('res-ppstar').textContent = fmt(ppstar);
      document.getElementById('res-rhorhostar').textContent = fmt(rhorhostar);
    }

    document.getElementById("computeIS").addEventListener("click", fillIsentropicTable);


    function computeForces() {
    // 1. Get input values
    const speed = parseFloat(document.getElementById('speed').value);
    const rho = parseFloat(document.getElementById('rho').value);
    const area = parseFloat(document.getElementById('area').value);
    const cd = parseFloat(document.getElementById('cd').value);
    const cl = parseFloat(document.getElementById('cl').value); // Get lift coefficient

    // Check for valid inputs

        // Check for valid inputs
    if (isNaN(speed) || isNaN(rho) || isNaN(area)) {
        // Handle invalid inputs
        document.getElementById('q').textContent = '—';
        document.getElementById('drag').textContent = '—';
        document.getElementById('lift').textContent = '—'; // Clear lift output
        return;
    }

    // 2. Perform calculations
    const dynamicPressure = 0.5 * rho * speed * speed;
    const dragForce = dynamicPressure * area * cd;
    const liftForce = dynamicPressure * area * cl; // Calculate lift force

    // 3. Display results
    document.getElementById('q').textContent = dynamicPressure.toFixed(2);
    document.getElementById('drag').textContent = dragForce.toFixed(2);
    document.getElementById('lift').textContent = liftForce.toFixed(2); // Display lift force

    if (isNaN(cd)) {
        // Handle invalid inputs
        document.getElementById('drag').textContent = '—';
        return;
    }
        // Check for valid inputs
    if (isNaN(cl)) {
        // Handle invalid inputs
        document.getElementById('lift').textContent = '—'; // Clear lift output
        return;
    }
}

// Add an event listener to call the function when the button is clicked
document.getElementById('computedrag').addEventListener('click', computeForces);