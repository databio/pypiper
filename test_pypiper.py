#!/usr/env python
"""
These are unit tests for the Pypiper class. Run with
python test_pypiper.py
"""

import unittest
import pypiper
import shutil #for rmtree
import os
import time
import subprocess

class PypiperTest(unittest.TestCase):

	def setUp(self):
		print("Setting up...")
		# Create a fixture
		self.pp = pypiper.Pypiper(name="sample_pipeline", outfolder="pipeline_output/", multi=True)
		self.pp2 = pypiper.Pypiper(name="sample_pipeline2", outfolder="pipeline_output/", multi=True)

	def tearDown(self):
		print("Tearing down...")
		self.pp.stop_pipeline()
		self.pp2.stop_pipeline()
		print("Removing " + self.pp.pipeline_outfolder)
		#shutil.rmtree(self.pp.pipeline_outfolder)
		#shutil.rmtree(self.pp2.pipeline_outfolder)
		del self.pp

	def test_me(self):
		print("Testing initialization...")
		self.assertEqual(self.pp.pipeline_name, "sample_pipeline")
		self.assertEqual(self.pp2.pipeline_name, "sample_pipeline2")
		# it creates an outfolder
		self.assertTrue(os.path.exists(self.pp.pipeline_outfolder))

		print("Testing status flags...")
		self.pp.set_status_flag("testing")
		self.assertTrue(os.path.isfile(self.pp.pipeline_outfolder + "sample_pipeline_testing.flag"))
		self.pp.set_status_flag("running")
		self.assertFalse(os.path.isfile(self.pp.pipeline_outfolder + "sample_pipeline_testing.flag"))
		self.assertTrue(os.path.isfile(self.pp.pipeline_outfolder + "sample_pipeline_running.flag"))

		print("Testing waiting for locks...")
		self.pp2.wait=False
		self.pp.wait=False
		sleep_lock = self.pp.pipeline_outfolder + "lock.sleep"
		subprocess.Popen("sleep 2; rm " + sleep_lock, shell=True)
		self.pp.create_file(sleep_lock)
		print("Putting lock file: " + sleep_lock)
		cmd = "echo hello"
		stamp = time.time()
		self.pp.call_lock(cmd, lock_name="sleep")
		print("Elapsed: " + str(self.pp.time_elapsed(stamp)))
		self.assertTrue(self.pp.time_elapsed(stamp) > 2)
		print("Wait for subprocess...")
		self.pp.wait_for_process(self.pp.running_subprocess)
		self.pp2.wait=True
		self.pp.wait=True

		print("Make sure the pipeline respects files already existing...")
		target = self.pp.pipeline_outfolder + "tgt"
		self.pp.call_lock("echo first > " + target, target, shell=True)
		self.pp.call_lock("echo second > " + target, target, shell=True) # Should not run
		with open(target) as f:
			lines = f.readlines()
		self.assertEqual(lines, ['first\n'])

		print("Execute a targetless command...")
		self.pp.call_lock("echo third > " + target, target=None, lock_name="test", shell=True) #
		with open(target) as f:
			lines = f.readlines()
		self.assertEqual(lines, ['third\n'])

		print("Test intermediate file cleanup...")
		tgt1 = self.pp.pipeline_outfolder + "tgt1.temp"
		tgt2 = self.pp.pipeline_outfolder + "tgt2.temp"
		self.pp.call_lock("touch " + tgt1, tgt1)
		self.pp.call_lock("touch " + tgt2, tgt2)
		self.pp.cleanup_append(tgt1)
		self.pp.cleanup()
		self.assertFalse(os.path.isfile(tgt1))
		self.assertTrue(os.path.isfile(tgt2))
		# How to test waiting for locks?

if __name__ == '__main__':
	unittest.main()
