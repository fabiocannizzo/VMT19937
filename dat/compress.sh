#!/bin/sh

matname=$1

7za a -t7z -m0=lzma2 -mx=5 -md=16m ${matname}.7z ${matname}.bits


