"""pytest configuration for the ACCeL v0.3 baseline suite.

The baseline suite pins the *current* observable behaviour of Box / System /
Systems so that v2.0 internals can be introduced without silently changing
results for existing users. Tests marked ``quirk`` pin behaviour that looks
questionable; change them only together with an explicit decision.
"""

from accel.util.log import Log

# ACCeL logs every mutation at DEBUG level to stderr by default; keep test output readable.
Log.console(show=False)
