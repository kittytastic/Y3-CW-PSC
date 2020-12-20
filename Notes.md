# vTune
## Compile
```
g++ -O3 -fopenmp -fopt-info-vec -g ../Solutions/base_af.cpp
```

```
icpc -O3 -fopenmp -qopt-report-phase=vec -qopt-report=2 -xhost --std=c++0x -g ./step2.cpp
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
