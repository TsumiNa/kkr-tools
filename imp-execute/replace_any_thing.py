#! /usr/bin/python3
# Filenamr : replace-any_thing.py
# -*-coding : utf-8-*-


where = ""
which = "input_relax-inward-Pd"
what  = r"7685949"
to    = r"7686369"


try:
    import os
    import re as regex
except ImportError:
    print('ImportError')


def replace():
    if where is "":
        _where = os.path.dirname(__file__)
        print('Search root is: ' + _where)
    elif os.path.exists('where'):
        _where = where
        print('Search root is: ' + _where)
    else:
        print('The path is not exist')
        return
    if which is "":
        print('Need to specify a file name')
        return
    if what is "":
        print('Need to specify what you want to find out')
        return
    if to is "":
        print('Need to specify what you want to be replaced')
        return
    for root, dirs, files in os.walk(_where):
        for _file_name in files:
            if _file_name == which:
                _file = os.path.join(root, _file_name)
                print(_file + 'will be modified')
                _fp = open(_file, 'r')
                _content = _fp.read()
                _fp.close()
                _content = regex.sub(what, to, _content)
                _fp = open(_file, 'w')
                _fp.write(_content)
                _fp.close


def main():
    replace()
    os.system('pause')


if __name__ == '__main__':
    main()
