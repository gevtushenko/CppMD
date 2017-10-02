# CppMD
This project is simplification of [SimpleMD](https://github.com/GiovanniBussi/simplemd) project - molecular dynamics software
which can be used to simulate Lennard-Jones systems. 

## How to use
Build CppMD
```bash
mkdir build && cd build
cmake .. && make -j
```

Run CppMD

```
./CppMD -c FILE.xyz
```

Visualize results with Paraview by using point gaussian

![cppmdparaview](https://user-images.githubusercontent.com/9890394/31095721-80a5a5ea-a7c2-11e7-8cb0-321ba45a4738.PNG)


## Tests
All tests using boost unit test framework. How to run the test:

```
make test
```

How to run the test coverage:
```
lcov --zerocounters --directory .
lcov --capture --initial --directory . --output-file app

make test

lcov --no-checksum --directory . --capture --output-file app.info
genhtml app.info
```

## Docs
How to build documentation:
```
make doc
```

