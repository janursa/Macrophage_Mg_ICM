#!/bin/bash
# git config --global http.proxy http[s]://nourisa:start.123@strand.fzg.local:1080
#git config --global http.proxy http[s]://nourisa:Hamburg_145@vpn.hzg.de
# git config --global http.proxy 127.0.0.1:1080
git add -A
git commit  --amend -m $1
git push --force 
