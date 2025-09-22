#!/usr/bin/env python
"""PipelineManager tests."""

import glob
import os
import shutil
import subprocess
import time
import unittest

import pypiper
from pypiper.exceptions import SubprocessError
from pypiper.utils import pipeline_filepath

__author__ = "Nathan Sheffield"
__email__ = "nathan@code.databio.org"


class PipelineManagerTests(unittest.TestCase):
    """Tests for pypiper's PipelineManager."""

    OUTFOLDER = "tests/Data/pipeline_output"

    @classmethod
    def _clean(cls):
        for d in glob.glob(cls.OUTFOLDER + "*"):
            if os.path.isdir(d):
                print("Removing " + d)
                shutil.rmtree(d)

    def setUp(self):
        """Start each test case with two pipeline managers."""
        print("Setting up...")
        # Create a fixture
        self.pp = pypiper.PipelineManager(
            "sample_pipeline", outfolder=self.OUTFOLDER, multi=True
        )
        self.pp2 = pypiper.PipelineManager(
            "sample_pipeline2", outfolder=self.OUTFOLDER, multi=True
        )

        self.pp3 = pypiper.PipelineManager(
            "sample_pipeline3", outfolder=self.OUTFOLDER + "3", multi=True
        )

    def tearDown(self):
        """Scrub the decks after each test case completes."""
        print("Tearing down...")
        self.pp.stop_pipeline()
        self.pp2.stop_pipeline()
        self.pp3.stop_pipeline()
        print("Removing " + self.pp.outfolder)
        # shutil.rmtree(self.pp.outfolder)
        # shutil.rmtree(self.pp3.outfolder)
        self._clean()
        del self.pp
        del self.pp2
        del self.pp3

    def _isFile(self, filename):
        """Determine if the first manager has this file."""
        filepath = pipeline_filepath(self.pp, filename=filename)
        return os.path.isfile(filepath)

    def _assertFile(self, filename):
        """Assert that the named file exists for first pipeline manager."""
        try:
            assert self._isFile(filename)
        except AssertionError:
            outfolder_contents = os.listdir(self.pp.outfolder)
            print(
                "Pipeline outfolder contents:\n{}".format("\n".join(outfolder_contents))
            )
            raise

    def _assertNotFile(self, filename):
        """Assert that given file doesn't exist for first manager."""
        assert not self._isFile(filename)

    def _assertLines(self, expected, observed):
        """Assert equality between collections of lines."""
        if isinstance(observed, str) and os.path.isfile(observed):
            with open(observed, "r") as f:
                observed = f.readlines()
        self.assertListEqual(expected, [l.strip() for l in observed])

    @classmethod
    def tearDownClass(cls):
        """Ensure folder/file cleanup upon test class completion."""
        cls._clean()

    def test_me(self):
        print("Testing initialization...")

        # Names
        self.assertEqual(self.pp.name, "sample_pipeline")
        self.assertEqual(self.pp2.name, "sample_pipeline2")

        # Outfolder creation
        self.assertTrue(os.path.isdir(self.pp.outfolder))

        print("Testing status flags...")
        self.pp._set_status_flag("completed")
        self._assertFile("sample_pipeline_DEFAULT_SAMPLE_NAME_completed.flag")
        self.pp._set_status_flag("running")
        self._assertNotFile("sample_pipeline_DEFAULT_SAMPLE_NAME_testing.flag")
        self._assertFile("sample_pipeline_DEFAULT_SAMPLE_NAME_running.flag")

        print("Testing waiting for locks...")
        self.pp2.wait = False
        self.pp.wait = False
        sleep_lock = pipeline_filepath(self.pp, filename="lock.sleep")
        subprocess.Popen("sleep .5; rm " + sleep_lock, shell=True)
        self.pp._create_file(sleep_lock)
        print("Putting lock file: " + sleep_lock)
        cmd = "echo hello"
        stamp = time.time()
        self.pp.run(cmd, lock_name="sleep")
        print("Elapsed: " + str(self.pp.time_elapsed(stamp)))
        self.assertTrue(self.pp.time_elapsed(stamp) > 1)

        print("Wait for subprocess...")
        for p in self.pp.running_procs.copy():
            self.pp._wait_for_process(self.pp.running_procs[p]["p"])
        self.pp2.wait = True
        self.pp.wait = True

        print("Make sure the pipeline respects files already existing...")
        target = pipeline_filepath(self.pp, filename="tgt")
        if os.path.isfile(target):  # for repeat runs.
            os.remove(target)

        self.pp.run("echo first > " + target, target, shell=True)
        # Should not run
        self.pp.run("echo second > " + target, target, shell=True)
        with open(target) as f:
            lines = f.readlines()
        self._assertLines(["first"], lines)

        print("Execute a targetless command...")
        self.pp.run("echo third > " + target, target=None, lock_name="test", shell=True)
        with open(target) as f:
            lines = f.readlines()
        self._assertLines(["third"], lines)

        # Test reporting results
        self.pp.report_result("key1", "abc")
        self.pp.report_result("key2", "def")
        key1 = self.pp.get_stat("key1")
        self.assertEqual(key1, "abc")

        try:
            key1 = self.pp2.get_stat("key1")  # should fail
        except KeyError:
            key1 = None
        self.assertEqual(key1, None)
        # We can no longer group based on 'shared' annotations.
        # key2 = self.pp2.get_stat("key2")  # should succeed
        # self.assertEqual(key2, "def")

        print("Test intermediate file cleanup...")
        tgt1 = pipeline_filepath(self.pp, filename="tgt1.temp")
        tgt2 = pipeline_filepath(self.pp, filename="tgt2.temp")
        tgt3 = pipeline_filepath(self.pp, filename="tgt3.temp")
        tgt4 = pipeline_filepath(self.pp, filename="tgt4.txt")
        tgt5 = pipeline_filepath(self.pp, filename="tgt5.txt")
        tgt6 = pipeline_filepath(self.pp, filename="tgt6.txt")
        tgt8 = pipeline_filepath(self.pp, filename="tgt8.cond")
        tgt9 = pipeline_filepath(self.pp, filename="tgt9.cond")
        tgt10 = pipeline_filepath(self.pp, filename="tgt10.txt")

        self.pp.run(
            "touch "
            + tgt1
            + " "
            + tgt2
            + " "
            + tgt3
            + " "
            + tgt4
            + " "
            + tgt5
            + " "
            + tgt6,
            lock_name="test",
        )
        self.pp.run("touch " + tgt8 + " " + tgt9, lock_name="test")

        # In global dirty mode, even non-manual clean files should not be deleted:
        self.pp.dirty = True
        self.pp.clean_add(tgt3)
        self.pp.clean_add(pipeline_filepath(self.pp, filename="*.temp"))
        self.pp.clean_add(tgt4)
        self.pp.clean_add(tgt5, conditional=True)
        self.pp.clean_add(
            pipeline_filepath(self.pp, filename="*.cond"), conditional=True
        )
        self.pp._cleanup()

        self.assertTrue(os.path.isfile(tgt1))
        self.assertTrue(os.path.isfile(tgt2))
        self.assertTrue(os.path.isfile(tgt3))
        self.assertTrue(os.path.isfile(tgt4))

        # But, these files *should* have been added to the cleanup script

        # test from different wd:
        tgt6_abs = os.path.abspath(tgt6)
        cfile = self.pp.cleanup_file
        ofolder = self.pp.outfolder
        cwd = os.getcwd()
        self.pp.clean_add(tgt6_abs)

        os.chdir("tests/Data/pipeline_output")
        self.pp.outfolder = "../" + ofolder
        self.pp.cleanup_file = "../" + cfile
        self.pp.clean_add(tgt6_abs)
        os.chdir(cwd)
        self.pp.cleanup_file = cfile
        self.pp.outfolder = ofolder

        print("Test manual cleanup adds")
        with open(self.pp.cleanup_file) as f:
            lines = f.readlines()

        print(lines)

        self.assertTrue(lines[2] == "rm tgt3.temp\n")
        self.assertTrue(lines[10] == "rm tgt6.txt\n")
        # lines is only 0-10 so the below code will error.
        # self.assertTrue(lines[11] == "rm tgt6.txt\n")

        self.pp.report_result("Test figure", os.path.join("fig", "fig.jpg"))

        # But in regular mode, they should be deleted:
        self.pp.dirty = False
        self.pp.clean_add(pipeline_filepath(self.pp, filename="*.temp"))
        self.pp.clean_add(tgt4)
        self.pp.clean_add(tgt5, conditional=True)
        self.pp.clean_add(
            pipeline_filepath(self.pp, filename="*.cond"), conditional=True
        )
        self.pp._cleanup()

        self.assertFalse(os.path.isfile(tgt1))
        self.assertFalse(os.path.isfile(tgt2))
        self.assertFalse(os.path.isfile(tgt3))
        self.assertFalse(os.path.isfile(tgt4))

        tgt7 = pipeline_filepath(self.pp, filename="tgt7.txt")
        self.pp.run("touch " + tgt7, tgt7)
        self.pp.clean_add(tgt7, manual=True)

        self.pp.run("touch " + tgt10, target=tgt10, clean=True)

        # Conditional delete should not delete tgt5
        # while pp2 is running
        self.assertTrue(os.path.isfile(tgt5))
        self.assertTrue(os.path.isfile(tgt8))
        self.assertTrue(os.path.isfile(tgt9))
        self.assertTrue(os.path.isfile(tgt10))  # auto cleanup

        # Stopping pp2 should cause tgt5 to be deleted
        self.pp2.stop_pipeline()
        self.pp._cleanup()
        self.assertFalse(os.path.isfile(tgt5))
        self.assertFalse(os.path.isfile(tgt8))
        self.assertFalse(os.path.isfile(tgt9))
        self.assertFalse(os.path.isfile(tgt10))

        # Manual clean should not clean
        self.assertTrue(os.path.isfile(tgt7))

        # cleanup should run on termination:
        self.pp.run("touch " + tgt6, tgt6)
        self.pp.clean_add(tgt6, conditional=True)
        self.pp.stop_pipeline()
        self.assertFalse(os.path.isfile(tgt5))

        # Manual clean should not clean even after pipeline stops.
        self.assertTrue(os.path.isfile(tgt7))

        print("Test failure and nofail options...")

        cmd = "thiscommandisbad"

        # Should not raise an error
        self.pp.run(cmd, target=None, lock_name="badcommand", nofail=True)
        self.pp.callprint(cmd, shell=None, nofail=True)

        # Should raise an error
        with self.assertRaises(SubprocessError):
            self.pp.run(cmd, target=None, lock_name="badcommand")

        print("Test dynamic recovery...")
        # send sigint
        self.pp.locks.append("lock.sleep")
        with self.assertRaises(KeyboardInterrupt):
            self.pp._signal_int_handler(None, None)

        sleep_lock = pipeline_filepath(self.pp, filename="lock.sleep")
        # subprocess.Popen("sleep .5; rm " + sleep_lock, shell=True)
        self.pp._create_file(sleep_lock)
        cmd = "echo hello"
        self.pp.run(cmd, lock_name="sleep")

        # subprocess.Popen("sleep .5; rm " + sleep_lock, shell=True)

        print("Test new start")
        if os.path.isfile(target):  # for repeat runs.
            os.remove(target)
        self.pp.run("echo first > " + target, target, shell=True)
        # Should not run
        self.pp.run("echo second > " + target, target, shell=True)
        with open(target) as f:
            lines = f.readlines()
        self._assertLines(["first"], lines)
        self.pp.new_start = True
        # Should run
        self.pp.run("echo third > " + target, target, shell=True)
        with open(target) as f:
            lines = f.readlines()
        self._assertLines(["third"], lines)

        print("Test dual target")
        self.pp.new_start = False
        if os.path.isfile(tgt1):
            os.remove(tgt1)
        self.pp.run("touch " + tgt6, tgt6)
        self.assertTrue(os.path.isfile(tgt6))
        # if target exists, should not run
        self.pp.run("touch " + tgt1, tgt6)
        self.assertFalse(os.path.isfile(tgt1))
        # if two targets, only one exists, should run
        self.assertFalse(os.path.isfile(tgt5))
        self.pp.run("touch " + tgt1, [tgt5, tgt6])
        self.assertTrue(os.path.isfile(tgt1))
        # if two targets, both exist, should not run
        self.assertFalse(os.path.isfile(tgt5))
        self.pp.run("touch " + tgt5, [tgt1, tgt6])
        self.assertFalse(os.path.isfile(tgt5))
        self.pp.pipestat.clear_status(self.pp.name, flag_names=["failed"])
        self.pp2.pipestat.clear_status(self.pp2.name, flag_names=["failed"])
        self.pp3.pipestat.clear_status(self.pp3.name, flag_names=["failed"])


def _make_pipe_filepath(pm, filename):
    return os.path.join(pm.outfolder, filename)


if __name__ == "__main__":
    unittest.main()
