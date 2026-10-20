"""Real-WASM regressions for method changes, presets, scans, late replies and reopening."""
from browser_controls import control

import argparse,functools,http.server,json,threading
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--browser-path');args=parser.parse_args()
class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self,*args):pass
    def end_headers(self):
        self.send_header('Cache-Control','no-store')
        super().end_headers()
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),functools.partial(Quiet,directory=str(ROOT/'build/wasm/dist')))
threading.Thread(target=server.serve_forever,daemon=True).start()
BRIDGE="""(() => {
  const Native=window.Worker;window.__requests=[];window.__workers=[];
  window.Worker=class {
    constructor(...args) {
      this.native=new Native(...args);window.__workers.push(this);
      this.native.onmessage=event=>{
        const data=window.__badBuild&&event.data.type==='ready'?{...event.data,buildId:'old-build'}:event.data;
        if(data.type==='result') {
          this.lastResult=data;
          setTimeout(()=>this.onmessage?.({data}),350);
        } else this.onmessage?.({data});
      };
      this.native.onerror=e=>this.onerror?.(e);
    }
    postMessage(data,...rest) {if(data.type!=='cancel')window.__requests.push(data);this.native.postMessage(data,...rest);}
    terminate(){this.native.terminate();}
  };
})();"""
try:
 with sync_playwright() as p:
    launch={'headless':True,'args':['--enable-unsafe-swiftshader']}
    if args.browser_path:launch['executable_path']=args.browser_path
    browser=p.chromium.launch(**launch);context=browser.new_context(viewport={'width':1512,'height':982},reduced_motion='reduce');context.add_init_script(BRIDGE)
    page=context.new_page();errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    url=f'http://127.0.0.1:{server.server_port}'
    def status(prefix):page.wait_for_function('(p)=>document.querySelector("#status").textContent.startsWith(p)',arg=prefix,polling=25)
    def run():page.keyboard.press('Control+Enter')
    def result():
        with page.expect_download() as d:control(page, '#download').click()
        return json.loads(Path(d.value.path()).read_text())
    page.goto(url);status('CONVERGED')
    # An explicit method must survive an example that recommends a different one.
    for method in ['mp2','uhf','rhf','mp2','uhf','mp2']:
        page.locator('#method').select_option(method);control(page, '#preset').select_option('hehplus')
        page.wait_for_timeout(180)
        assert page.locator('#method').input_value()==method
        run();status('CONVERGED');r=result();assert r['method']==method.upper()
        assert page.evaluate('window.__requests.at(-1).options.method')==method
        if method=='mp2':assert r['mp2']['status']=='complete'
    # The stored molecule/method/charge/spin must restore together on reload.
    page.reload();status('CONVERGED')
    assert page.locator('#method').input_value()=='mp2'
    assert control(page, '#charge').input_value()=='1'
    assert control(page, '#preset').input_value()=='hehplus'
    assert result()['method']=='MP2'
    # Reopening another tab at the same origin restores the same coherent draft.
    reopened=context.new_page();reopened.goto(url)
    reopened.wait_for_function('document.querySelector("#status").textContent.startsWith("CONVERGED")')
    assert reopened.locator('#method').input_value()=='mp2'
    assert reopened.locator('#charge').input_value()=='1';reopened.close()
    # Keep a finished RHF scan, then start MP2. The delayed real result provides
    # a deterministic window in which clicking the old scan used to reset it.
    page.locator('#method').select_option('rhf');control(page, '#preset').select_option('h2');control(page, '#mode').select_option('scan')
    for key,value in [('a','1'),('b','2'),('moving','2'),('start','.7'),('end','1.3'),('count','3')]:control(page, '#scan-'+key).fill(value)
    run();page.wait_for_function('document.querySelector("#progress-text").textContent.startsWith("Complete")')
    control(page, '#mode').select_option('rhf');page.locator('#method').select_option('mp2')
    page.locator('#scan-point-table tbody tr').first.click()
    assert page.locator('#method').input_value()=='mp2'
    assert 'current setup MP2/' in page.locator('#status').inner_text()
    run()
    assert page.locator('#method').is_disabled()
    page.locator('#tab-pes').click();page.locator('#scan-point-table tbody tr').first.click()
    assert page.locator('#method').input_value()=='mp2'
    status('CONVERGED');assert result()['method']=='MP2'
    # A completed result delivered again after a settings edit must be ignored.
    page.locator('#method').select_option('uhf')
    page.evaluate('window.__workers[0].onmessage({data:window.__workers[0].lastResult})')
    assert 'Settings changed' in page.locator('#status').inner_text()
    assert control(page, '#download').is_disabled()
    # Test native keyboard selection after the former debounce interval.
    page.locator('#method').focus();page.keyboard.press('Home');page.wait_for_timeout(180);page.keyboard.press('End');page.keyboard.press('Enter')
    assert page.locator('#method').input_value()=='mp2'
    run();status('CONVERGED');assert result()['method']=='MP2'
    bad=context.new_page();bad.add_init_script('window.__badBuild=true');bad.goto(url)
    bad.wait_for_function('document.querySelector("#status").textContent.includes("versions differ")')
    assert bad.locator('#calculate').is_disabled();bad.close()
    assert not errors,errors
    browser.close();print('Repeated method/preset changes, active-job scan guard, stale replies, reload/reopen, native selection and build mismatch PASS')
finally:server.shutdown();server.server_close()
