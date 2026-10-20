"""Check complete module-graph versioning without requiring an Emscripten build."""
import importlib.util,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('package_web',ROOT/'scripts/package_web.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class PackageTests(unittest.TestCase):
    def test_consistent_rebuild_and_changed_assets(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);web=root/'web';dist=root/'dist';web.mkdir();dist.mkdir()
            (web/'index.html').write_text('<link href="./styles.css"><script src="./vendor/viewer.js"></script><script type="module" src="./app.mjs"></script>')
            (web/'app.mjs').write_text("import './settings.mjs'; new Worker(new URL('./worker.mjs',import.meta.url),{type:'module'});")
            (web/'worker.mjs').write_text("import './build-info.mjs'; import './qc_core.mjs';")
            (web/'settings.mjs').write_text('export const schema=2;')
            (web/'build-info.mjs').write_text("export const BUILD_ID='source';")
            (web/'styles.css').write_text('body {color:green}')
            (web/'vendor').mkdir();(web/'vendor/viewer.js').write_text('window.viewer=true;')
            (dist/'qc_core.mjs').write_text('export default ()=>{};');(dist/'qc_core.wasm').write_bytes(b'wasm-fixture')
            first=module.package(web,dist)
            for path,ref in [('index.html','./app.mjs'),('index.html','./vendor/viewer.js'),('app.mjs','./worker.mjs'),('app.mjs','./settings.mjs'),('worker.mjs','./qc_core.mjs')]:
                self.assertIn(ref+'?v='+first,(dist/path).read_text())
            self.assertIn("BUILD_ID='"+first+"'",(dist/'build-info.mjs').read_text())
            self.assertEqual(module.package(web,dist),first)
            (web/'styles.css').write_text('body {color:blue}')
            second=module.package(web,dist);self.assertNotEqual(second,first)
            self.assertIn('./worker.mjs?v='+second,(dist/'app.mjs').read_text())
            (dist/'qc_core.wasm').write_bytes(b'new-wasm-fixture')
            self.assertNotEqual(module.package(web,dist),second)

if __name__=='__main__':unittest.main()
