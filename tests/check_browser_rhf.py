"""Optional RHF browser test: python3 tests/check_browser_rhf.py [--browser-path PATH]."""
import argparse
import functools
import http.server
import json
from pathlib import Path
import tempfile
import threading
import numpy as np
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--browser-path')
args=parser.parse_args()
reference=json.loads((ROOT/'tests/fixtures/water-rhf-sto3g.json').read_text())
handler=functools.partial(http.server.SimpleHTTPRequestHandler,directory=str(ROOT/'build/wasm/dist'))
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),handler)
threading.Thread(target=server.serve_forever,daemon=True).start()
try:
    with sync_playwright() as p,tempfile.TemporaryDirectory() as temp:
        options={'headless':True}
        if args.browser_path: options['executable_path']=args.browser_path
        browser=p.chromium.launch(**options)
        page=browser.new_page(viewport={'width':1280,'height':1000})
        errors=[]
        page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(f'http://127.0.0.1:{server.server_port}')
        def wait_status(prefix):
            page.wait_for_function('(prefix) => document.querySelector("#status").textContent.startsWith(prefix)',arg=prefix,timeout=30000)
        def download():
            with page.expect_download() as event: page.locator('#download').click()
            file=Path(temp)/'rhf.json'; event.value.save_as(file)
            return json.loads(file.read_text())
        wait_status('CONVERGED')
        assert '-74.963063129729' in page.locator('#summary').inner_text()
        data=download(); actual=data['rhf']
        assert actual['converged']
        assert abs(actual['energy']-reference['energy'])<1e-9
        for key in ['density','fock','orbitalEnergies']:
            assert np.max(np.abs(np.asarray(actual[key])-reference[key]))<1e-7,key
        assert page.locator('#scf-history tbody tr').count()==actual['iterations']
        assert page.locator('#orbital-table tr').count()==8
        for kind in ['density','fock']:
            page.locator('#matrix-kind').select_option(kind)
            assert page.locator('#matrix td').count()==49
        # Observe progress during a fresh run, not just the final history render.
        page.evaluate('''() => {
          window.iterationsObserved=0;
          new MutationObserver(() => {if(document.querySelector('#status').textContent.startsWith('SCF iteration')) window.iterationsObserved++;})
            .observe(document.querySelector('#status'),{childList:true});
        }''')
        page.locator('#calculate').click(); wait_status('CONVERGED')
        assert page.evaluate('window.iterationsObserved')>0
        page.locator('#maxiter').fill('1'); page.locator('#calculate').click()
        wait_status('NOT CONVERGED')
        limited=download()['rhf']
        assert not limited['converged'] and limited['status']=='not_converged'
        assert 'Unconverged energy' in page.locator('#summary').inner_text()
        page.locator('#maxiter').fill('100'); page.locator('#charge').fill('1')
        assert page.locator('#download').is_disabled()
        page.locator('#calculate').click(); wait_status('Input or calculation error')
        assert 'even electron count' in page.locator('#summary').inner_text()
        page.locator('#water').click(); wait_status('CONVERGED')
        page.locator('#geometry').fill('H 0 0 0\nH 0 0 1.4')
        page.locator('#unit').select_option('bohr'); page.locator('#calculate').click()
        wait_status('CONVERGED')
        assert '-1.116714325063' in page.locator('#summary').inner_text()
        page.locator('#water').click(); wait_status('CONVERGED')
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth<=window.innerWidth')
        assert not errors,errors
        browser.close()
        print('RHF browser energy, density, orbitals, JSON, live iterations: PASS')
        print('Nonconvergence, charge rejection, recovery, H2/bohr, mobile layout: PASS')
finally:
    server.shutdown(); server.server_close()
