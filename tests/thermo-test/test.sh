#!/bin/bash

BUILD_DIR=$1

$BUILD_DIR/apollo -C config/config.toml -S tests/thermo-test/thermotest.toml
