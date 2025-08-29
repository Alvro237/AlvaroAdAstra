    // Utility: format numbers nicely
    const fmt = (x) => Number.isFinite(x) ? x.toLocaleString(undefined, {maximumFractionDigits: 3}) : '—';

    function compute(){
      const V   = parseFloat(document.querySelector('#speed').value);
      const rho = parseFloat(document.querySelector('#rho').value);
      const A   = parseFloat(document.querySelector('#area').value);
      const Cd  = parseFloat(document.querySelector('#cd').value);

      // q = 0.5 * rho * V^2
      const q = 0.5 * rho * V * V;
      // D = q * Cd * A
      const D = q * Cd * A;

      document.querySelector('#q').textContent = fmt(q);
      document.querySelector('#drag').textContent = fmt(D);
      // Persist last values so they survive refreshes
      localStorage.setItem('calc-state', JSON.stringify({V, rho, A, Cd}));
    }

    // Hook up events
    document.querySelector('#compute').addEventListener('click', compute);
    ['#speed','#rho','#area','#cd'].forEach(sel => document.querySelector(sel).addEventListener('input', compute));

    // Restore last state if present
    try{
      const saved = JSON.parse(localStorage.getItem('calc-state'));
      if(saved){
        document.querySelector('#speed').value = saved.V;
        document.querySelector('#rho').value   = saved.rho;
        document.querySelector('#area').value  = saved.A;
        document.querySelector('#cd').value    = saved.Cd;
      }
    }catch{}

    // Set current year and do an initial compute
    document.querySelector('#year').textContent = new Date().getFullYear();
    compute();
