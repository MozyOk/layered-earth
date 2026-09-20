import unittest
import tempfile
import subprocess
import sys
import json
from pathlib import Path
import numpy as np
from layered_earth import analyze, R_KM

class ModelTests(unittest.TestCase):
    def run_grid(self,h,**kw):
        h=np.array(h,dtype=float)
        return analyze(h,np.arange(h.shape[1]+1),np.arange(h.shape[0]+1),**kw)

    def test_half_open_negative_and_missing(self):
        layers,_,c,_,_=self.run_grid([[-.1,0,4.999,5,np.nan]])
        np.testing.assert_equal(layers,[[-1,0,0,1,np.nan]])
        self.assertEqual(sum(x['cells'] for x in c),4)

    def test_diagonal(self):
        h=[[1,6],[6,1]]
        self.assertEqual(len(self.run_grid(h,connectivity=4)[2]),4)
        self.assertEqual(len(self.run_grid(h,connectivity=8)[2]),2)

    def test_area_conservation(self):
        _,_,c,_,_=self.run_grid([[0,0],[5,10]])
        expected=R_KM**2*np.deg2rad(2)*np.sin(np.deg2rad(2))
        self.assertAlmostEqual(sum(x['area_km2'] for x in c),expected)

    def test_upward_contacts_allow_skipped_layers(self):
        _,_,_,_,e=self.run_grid([[10,0,5]])
        self.assertEqual(sorted(x['layer_jump'] for x in e),[1,2])
        self.assertEqual(len(e),2)

    def test_mask(self):
        self.assertEqual(len(self.run_grid([[0,0,0]],mask=[[1,0,1]])[2]),2)

    def test_empty(self):
        with self.assertRaises(ValueError): self.run_grid([[np.nan]])

    def test_invalid_band(self):
        for band in [0,-1,np.nan]:
            with self.assertRaises(ValueError): self.run_grid([[0]],band=band)

    def test_global_rejected(self):
        with self.assertRaises(ValueError): analyze([[0]],[-180,180],[-90,90])

    def test_cli_exports(self):
        # Synthetic fixture verifies rendering only; never published as Earth data.
        with tempfile.TemporaryDirectory() as root:
            root=Path(root)
            np.savez(root/'fixture.npz',h=[[0.,5.],[10.,15.]],lon_edges=[139,140,141],lat_edges=[35,36,37])
            subprocess.run([sys.executable,'layered_earth.py','--dem',str(root/'fixture.npz'),
                            '--source','synthetic-unit-test-only','--output',str(root/'out')],
                           check=True,capture_output=True)
            meta=json.loads((root/'out/provenance.json').read_text())
            self.assertEqual(meta['components'],4)
            for name in ['overview.svg','overview.png','results.md','components.csv','contacts.csv','layers.csv']:
                self.assertGreater((root/'out'/name).stat().st_size,0)

if __name__=='__main__': unittest.main()
