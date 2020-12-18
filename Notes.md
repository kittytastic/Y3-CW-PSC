# vTune
## Compile
```
gcc -O3 -g ../Source/base_af.cpp
```


## Source
```
source /opt/intel/oneapi/setvars.sh
```

## Start server
```
 vtune-backend --web-port=8080 --allow-remote-ui
```

## Start agent
```
~/intel/ProfileAgent/vtune-agent --owner passphrase_authenticated_user
```
