#! /usr/bin/env python

"""This module creates a Timer class that can be used for timing pieces of code.
Inspired by a tutorial on RealPython: https://realpython.com/python-timer/
"""
import time


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""


class Timer:
    timers = dict()

    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self, name=None):
        """Stop the timer, and save the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        # Add new named timers to dictionary of timers
        if name:
            self.timers[name] = elapsed_time

        return elapsed_time

    def sum(self, key_str=''):
        """Sum the times for all dictionary entries where the key
        contains key_str"""
        total_time = 0.

        for key in self.timers:
            if key_str in key:
                total_time += self.timers[key]
        return total_time
