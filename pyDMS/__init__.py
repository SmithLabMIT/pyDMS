r'''
pyDMS.dms
                 ____  __  ________
    ____  __  __/ __ \/  |/  / ___/
   / __ \/ / / / / / / /|_/ /\__ \
  / /_/ / /_/ / /_/ / /  / /___/ /
 / .___/\__, /_____/_/  /_//____/
/_/    /____/

Copyright 2025 Brandon C. Tapia

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the “Software”),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import warnings

__version__ = '0.0.1dev1'


def error_in_red(message):
    '''
    Red error message
    '''
    print("\033[91m" + message + "\033[0m")
    warnings.warn(message, UserWarning)


def warning_in_orange(message):
    '''
    Orange warning message
    '''
    print("\033[38;2;255;165;0m" + message + "\033[0m")
    warnings.warn(message, UserWarning)
