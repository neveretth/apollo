#!/usr/bin/env bash

echo "Running Valgrind with full leak check and origin tracing..."

valgrind --leak-check=full --track-origins=yes ../build/apollo -S ../simulation/valgrind.toml -C ../config/config.toml 2>&1 | tee valgrind_output.log

if grep -q "ERROR SUMMARY: 0 errors" valgrind_output.log; then
    echo ":) Valgrind test passed: No memory errors detected."
    exit 0
else
    echo ":( Valgrind test failed: Memory issues detected. See valgrind_output.log for details."
    exit 1
fi
