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
		#self.pp.stop_pipeline()
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
		subprocess.Popen("sleep .5; rm " + sleep_lock, shell=True)
		self.pp.create_file(sleep_lock)
		print("Putting lock file: " + sleep_lock)
		cmd = "echo hello"
		stamp = time.time()
		self.pp.call_lock(cmd, lock_name="sleep")
		print("Elapsed: " + str(self.pp.time_elapsed(stamp)))
		self.assertTrue(self.pp.time_elapsed(stamp) > 1)
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
		tgt3 = self.pp.pipeline_outfolder + "tgt3.temp"
		tgt4 = self.pp.pipeline_outfolder + "tgt4.txt"
		tgt5 = self.pp.pipeline_outfolder + "tgt5.txt"
		tgt6 = self.pp.pipeline_outfolder + "tgt6.txt"
		tgt8 = self.pp.pipeline_outfolder + "tgt8.cond"
		tgt9 = self.pp.pipeline_outfolder + "tgt9.cond"

		self.pp.call_lock("touch " + tgt1 + " " + tgt2 + " " + tgt3 + " " + tgt4 + " " + tgt5, lock_name="test")
		self.pp.call_lock("touch " + tgt8 + " " + tgt9, lock_name="test")

		self.pp.clean_add(self.pp.pipeline_outfolder + "*.temp")
		self.pp.clean_add(tgt4)
		self.pp.clean_add(tgt5, conditional=True)
		self.pp.clean_add(self.pp.pipeline_outfolder +"*.cond", conditional=True)
		self.pp.cleanup()

		self.assertFalse(os.path.isfile(tgt1))
		self.assertFalse(os.path.isfile(tgt2))
		self.assertFalse(os.path.isfile(tgt3))
		self.assertFalse(os.path.isfile(tgt4))

		tgt7 = self.pp.pipeline_outfolder + "tgt7.txt"
		self.pp.call_lock("touch " + tgt7, tgt7)
		self.pp.clean_add(tgt7, manual=True)



		# Conditional delete should not delete tgt5
		# while pp2 is running
		self.assertTrue(os.path.isfile(tgt5))
		self.assertTrue(os.path.isfile(tgt8))
		self.assertTrue(os.path.isfile(tgt9))

		# Stopping pp2 should cause tgt5 to be deleted
		self.pp2.stop_pipeline()
		self.pp.cleanup()
		self.assertFalse(os.path.isfile(tgt5))
		self.assertFalse(os.path.isfile(tgt8))
		self.assertFalse(os.path.isfile(tgt9))

		# Manual clean should not clean
		self.assertTrue(os.path.isfile(tgt7))

		# cleanup should run on termination:
		self.pp.call_lock("touch " + tgt6, tgt6)
		self.pp.clean_add(tgt6, conditional=True)
		self.pp.stop_pipeline()
		self.assertFalse(os.path.isfile(tgt5))

		# Manual clean should not clean even after pipeline stops
		self.assertTrue(os.path.isfile(tgt7))




if __name__ == '__main__':
	unittest.main()
