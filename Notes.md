# vTune
## Compile
```
g++ -O3 -fopenmp -g ../Solutions/base_af.cpp
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
 wget https://172.29.146.40:8080/api/collection-agent/download --no-check-certificate
~/intel/ProfilerAgent/vtune-agent --owner passphrase_authenticated_user
```
