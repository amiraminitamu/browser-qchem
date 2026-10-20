"""Real-WASM integration checks for the canvas-first two-zone workbench."""
import argparse, functools, http.server, json, threading
from browser_controls import control
from pathlib import Path
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--browser-path')
args = parser.parse_args()

class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args): pass

server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(Quiet, directory=str(ROOT/'build/wasm/dist')))
threading.Thread(target=server.serve_forever, daemon=True).start()
try:
    with sync_playwright() as p:
        launch = {'headless': True, 'args': ['--enable-unsafe-swiftshader']}
        if args.browser_path: launch['executable_path'] = args.browser_path
        browser = p.chromium.launch(**launch)
        page = browser.new_page(viewport={'width': 1512, 'height': 900}, reduced_motion='reduce')
        errors = []
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.goto(f'http://127.0.0.1:{server.server_port}')

        def status(prefix):
            page.wait_for_function('(p)=>document.querySelector("#status").textContent.startsWith(p)', arg=prefix)

        def open_panel(id):
            panel = page.locator('#'+id)
            if not panel.evaluate('(e)=>e.open'): panel.locator(':scope > summary').click()

        def export_result():
            open_panel('export-popover')
            with page.expect_download() as event: page.locator('#download').click()
            return json.loads(Path(event.value.path()).read_text())

        def in_viewport():
            assert page.evaluate('document.documentElement.scrollWidth<=innerWidth && document.documentElement.scrollHeight<=innerHeight')
            for selector in ['#calculate', '.viewer-card', '#results-section']:
                box = page.locator(selector).bounding_box()
                assert box and box['x'] >= 0 and box['y'] >= 0, (selector, box)
                assert box['x']+box['width'] <= page.viewport_size['width']+1, (selector, box)
                assert box['y']+box['height'] <= page.viewport_size['height']+1, (selector, box)

        status('CONVERGED')
        assert page.locator('.sidebar,.activity-rail,.menu-bar,.app-toolbar').count() == 0
        assert page.locator('button[type=submit]').count() == 1
        assert page.locator('#job-count').inner_text() == '1'
        assert page.locator('#coordinates-editor').evaluate('(e)=>!e.open')
        original_energy = page.locator('#energy-value').inner_text()
        left = page.locator('.viewer-card').bounding_box()
        right = page.locator('#results-section').bounding_box()
        assert right['x'] >= left['x']+left['width']
        in_viewport()

        # Pointer and keyboard split resizing, including explicit bounds.
        handle = page.locator('#split-handle')
        box = handle.bounding_box()
        page.mouse.move(box['x']+box['width']/2, box['y']+box['height']/2)
        page.mouse.down(); page.mouse.move(640, 450); page.mouse.up()
        assert 40 <= int(handle.get_attribute('aria-valuenow')) <= 45
        handle.focus(); page.keyboard.press('Home')
        assert handle.get_attribute('aria-valuenow') == '35'
        page.keyboard.press('End'); assert handle.get_attribute('aria-valuenow') == '70'
        handle.dblclick(); assert handle.get_attribute('aria-valuenow') == '55'

        # Analysis controls and orbital selection preserve the numerical result.
        page.locator('#tab-orbitals').click()
        page.locator('#orbital-table tbody tr').first.click()
        assert page.locator('#orbital-overlay-slot').is_visible()
        assert page.locator('#orbital-table tr.selected').count() == 1
        page.locator('#tab-matrices').click(); page.locator('#matrix-kind').select_option('core')
        assert page.locator('#matrix tbody tr').count() == 7
        assert page.locator('#energy-value').inner_text() == original_energy
        assert export_result()['rhf']['converged']
        page.keyboard.press('Escape')
        page.locator('#tab-data').click()
        assert page.locator('#coordinate-table tbody tr').count() == 3
        assert page.locator('#distance-matrix tbody tr').count() == 3
        page.locator('#tab-scf').click(); page.keyboard.press('ArrowRight')
        assert page.locator('#tab-orbitals').get_attribute('aria-selected') == 'true'

        page.locator('#xyz-toggle').click()
        assert page.locator('#geometry').is_visible()
        page.keyboard.press('Escape')
        assert page.locator('#xyz-toggle').get_attribute('aria-expanded') == 'false'
        page.locator('#jobs-toggle').click(); assert page.locator('#job-list .job-item').count() == 1
        page.keyboard.press('Escape'); assert page.locator('#jobs-drawer').is_hidden()
        page.locator('#help-open').click(); assert page.locator('#workspace-guide').evaluate('(e)=>e.open')
        page.keyboard.press('Escape')
        page.locator('#focus-view').click(); assert page.locator('#results-section').is_hidden()
        page.keyboard.press('Escape'); assert page.locator('#results-section').is_visible()

        open_panel('solver-popover'); page.locator('#maxiter').fill('1')
        page.keyboard.press('Control+Enter'); status('NOT CONVERGED')
        assert page.locator('.summary-details').evaluate('(e)=>e.open')
        open_panel('solver-popover'); page.locator('#maxiter').fill('100')
        page.keyboard.press('Control+Enter'); status('CONVERGED')
        page.locator('#tab-scf').click()
        page.locator('.summary-details').evaluate('(e)=>e.open=false')
        page.keyboard.press('Escape')
        page.screenshot(path=str(ROOT.parent/'redesign-desktop.png'))

        # The same controls remain canonical on a phone.
        page.set_viewport_size({'width': 390, 'height': 844})
        open_panel('molecule-popover'); page.locator('#preset').select_option('h2')
        page.locator('#basis').select_option('6-31g')
        page.locator('#calculate').click(); status('CONVERGED')
        open_panel('task-popover'); page.locator('#mode').select_option('scan')
        for key, value in [('a','1'), ('b','2'), ('moving','2'), ('start','.5'), ('end','1.5'), ('count','3')]:
            control(page, '#scan-'+key).fill(value)
        page.locator('#calculate').click()
        page.wait_for_function('document.querySelector("#progress-text").textContent.startsWith("Complete")')
        page.locator('#scan-point-table tbody tr').nth(1).click()
        assert abs(float(page.locator('#energy-value').inner_text().replace('−', '-'))-(-1.094807962861)) < 1e-9
        energy = page.locator('#energy-value').inner_text()
        page.locator('#pes-energy-unit').select_option('eh')
        page.locator('#pes-relative').uncheck()
        assert page.locator('#energy-value').inner_text() == energy
        open_panel('export-popover')
        assert not page.locator('#scan-json').is_disabled()
        page.keyboard.press('Escape')
        in_viewport()
        page.screenshot(path=str(ROOT.parent/'redesign-mobile.png'))
        for width in [360, 390, 768, 801, 1024, 1280, 1512]:
            page.set_viewport_size({'width': width, 'height': 900}); in_viewport()
        assert not errors, errors
        browser.close()
        print('Canvas split, canonical controls, orbital selection, retained results, exports, drawers, solver failure/recovery and mobile PES PASS')
finally:
    server.shutdown(); server.server_close()
