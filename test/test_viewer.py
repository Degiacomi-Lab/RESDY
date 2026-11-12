import unittest
import sys
import os
from unittest.mock import patch
import dash

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD


EXPECTED_COMPONENT_IDS = {'tabs_analysis', 'tabs_content_analysis'}

class Test_Viewer(unittest.TestCase):
    def setUp(self):
        self.outdir = 'demo'
        self.df_input = f'{self.outdir}{os.sep}demo_input.csv'
        self.df_measures = f'{self.outdir}{os.sep}measures.csv'
        self.V = RD.Viewer(outdir=self.outdir,
                           df_measures=self.df_measures)
        self.run_name = 'run' if hasattr(self.V.data_viewer_app, 'run') else 'run_server'

    def test_app_is_dash_instance(self):
        self.assertIsInstance(self.V.data_viewer_app, dash.Dash)

    def test_layout_is_set(self):
        self.assertIsNotNone(self.V.data_viewer_app.layout)

    def test_layout_contains_expected_ids(self):
        layout = self.V.data_viewer_app.layout
        if callable(layout):
            layout = layout()
        ids = set(layout._traverse_ids()) if hasattr(layout, "_traverse_ids") else set()
        ids = {c.id for c in layout._traverse() if getattr(c, "id", None)} | ids
        for expected in EXPECTED_COMPONENT_IDS:
            self.assertIn(expected, ids)

    def test_callbacks_registered(self):
        self.assertGreater(len(self.V.data_viewer_app.callback_map), 0)

    def test_launch_viewer_calls_run_with_defaults(self):
        with patch.object(self.V.data_viewer_app, self.run_name) as mock_run:
            self.V.launch_viewer()
        mock_run.assert_called_once()

    def test_server_responds(self):
        '''Hit the Flask test client without actually launching the server.'''
        client = self.V.data_viewer_app.server.test_client()
        response = client.get("/")
        self.assertEqual(response.status_code, 200)


if __name__ == "__main__":
    unittest.main()
