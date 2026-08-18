r"""
pyDMS.dms
                 ____  __  ________
    ____  __  __/ __ \/  |/  / ___/
   / __ \/ / / / / / / /|_/ /\__ \
  / /_/ / /_/ / /_/ / /  / /___/ /
 / .___/\__, /_____/_/  /_//____/
/_/    /____/

Copyright 2026 Massachusetts Institute of Technology

Licensed under the BSD 3-Clause License

Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions
   and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of
   conditions and the following disclaimer in the documentation and/or other materials provided with
   the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to
   endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import warnings

# TODO: streamline versioning
__version__ = "0.9.0"


def error_in_red(message):
    """
    Red error message
    """
    print("\033[91m" + message + "\033[0m")
    warnings.warn(message, UserWarning)


def warning_in_orange(message):
    """
    Orange warning message
    """
    print("\033[38;2;255;165;0m" + message + "\033[0m")
    warnings.warn(message, UserWarning)
