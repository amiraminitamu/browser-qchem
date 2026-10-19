"""Optional end-to-end test: python3 tests/check_browser.py [--browser-path PATH].

Requires Playwright and its Chromium browser. The server is temporary and local.
"""
import argparse
import functools
import http.server
import json
from pathlib import Path
import tempfile
import threading

import numpy as np
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--browser-path')
args = parser.parse_args()
handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(ROOT/'build/wasm/dist'))
server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), handler)
threading.Thread(target=server.serve_forever, daemon=True).start()
reference = json.loads((ROOT/'tests/fixtures/water-sto3g.json').read_text())

try:
    with sync_playwright() as p, tempfile.TemporaryDirectory() as temp:
        options = {'headless':True,'args':['--enable-unsafe-swiftshader']}
        if args.browser_path:
            options['executable_path'] = args.browser_path
        browser = p.chromium.launch(**options)
        page = browser.new_page(viewport={'width':1280,'height':1000})
        errors = []
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.goto(f'http://127.0.0.1:{server.server_port}')
        page.wait_for_function("document.querySelector('#status').textContent.startsWith('CONVERGED')",timeout=30000)
        page.locator('#mode').select_option('integrals')
        page.locator('#calculate').click()
        def wait_pass():
            page.wait_for_function("document.querySelector('#status').textContent.startsWith('PASS')", timeout=30000)
        wait_pass()
        assert '2401 values' in page.locator('#summary').inner_text()
        assert page.locator('#matrix tr').count() == 8
        for kind in ['kinetic','nuclear','core','overlap']:
            page.locator('#matrix-kind').select_option(kind)
            assert page.locator('#matrix td:not(:first-child)').count() == 49
        with page.expect_download() as event:
            page.locator('#download').click()
        target = Path(temp)/'integrals.json'
        event.value.save_as(target)
        actual = json.loads(target.read_text())
        for key in ['overlap','kinetic','nuclear','eri']:
            assert np.max(np.abs(np.asarray(actual[key])-reference[key])) < 1e-10, key
        print('Browser worker and JSON export: all water integrals match PySCF PASS')
        page.locator('.eri-details summary').click()
        page.get_by_label('ERI index i',exact=True).fill('3')
        page.get_by_label('ERI index j',exact=True).fill('5')
        page.get_by_label('ERI index k',exact=True).fill('6')
        page.get_by_label('ERI index l',exact=True).fill('7')
        expected = reference['eri'][((2*7+4)*7+5)*7+6]
        assert f'{expected:.14f}' in page.locator('#eri-value').inner_text()
        page.get_by_label('ERI index i',exact=True).fill('0')
        assert 'valid integer' in page.locator('#eri-value').inner_text()
        page.locator('#geometry').fill('H 0 0 0\nH 0 0 1.4')
        assert page.locator('#matrix tr').count() == 0
        assert page.locator('#download').is_disabled()
        page.locator('#unit').select_option('bohr')
        page.locator('#calculate').click(); wait_pass()
        assert page.locator('#matrix tr').count() == 3
        assert '0.65931821' in page.locator('#matrix').inner_text()
        page.locator('#geometry').fill('H 0 0 0\nH 0 0 0')
        page.locator('#calculate').click()
        page.wait_for_function("document.querySelector('#status').textContent.includes('error')")
        assert 'coincide' in page.locator('#summary').inner_text()
        assert page.locator('#download').is_disabled()
        page.locator('#water').click(); wait_pass()
        assert page.locator('#matrix tr').count() == 8
        # Check narrow layout without horizontal page overflow.
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth')
        assert not errors, errors
        browser.close()
        print('Matrix/ERI inspection, units, invalid input, recovery, mobile layout: PASS')
finally:
    server.shutdown()
    server.server_close()
