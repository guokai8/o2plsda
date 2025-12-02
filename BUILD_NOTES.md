# Build Notes for o2plsda Package

## Architecture Compatibility Issue - RESOLVED

### Problem
The package was previously compiled for x86_64 (Intel) architecture, but you're running on ARM64 (Apple Silicon M-series chip). This caused the error:
```
unable to load shared object...
(mach-o file, but is an incompatible architecture (have 'x86_64', need 'arm64'))
```

### Solution Applied
1. **Removed old compiled files:**
   ```bash
   cd src/
   rm -f *.o *.so symbols.rds
   ```

2. **Rebuilt package for ARM64:**
   ```bash
   R CMD INSTALL --preclean --no-multiarch .
   ```

3. **Verification:**
   - Package now compiled as: `Mach-O 64-bit dynamically linked shared library arm64`
   - Package loads successfully: ✓
   - Core functions tested: ✓

### Current Status
**✓ FULLY FUNCTIONAL on ARM64 (Apple Silicon)**

---

## Building the Package

### For Development/Testing

If you make changes to the C++ code (*.cpp files in src/), you need to rebuild:

```bash
# Clean and rebuild
cd /path/to/o2plsda-main
rm -f src/*.o src/*.so src/symbols.rds
R CMD INSTALL --preclean .
```

### For R CMD check

```bash
# Build package tarball
R CMD build .

# Check package
R CMD check --as-cran o2plsda_0.0.28.tar.gz
```

### For CRAN Submission

```bash
# Clean everything
rm -rf src/*.o src/*.so src/symbols.rds

# Build source package
R CMD build .

# Final check
R CMD check --as-cran o2plsda_0.0.28.tar.gz
```

---

## Important Notes

### Architecture-Specific Compiled Files

The following files in `src/` are architecture-specific and should **NOT** be committed to git:
- `*.o` - Object files
- `*.so` - Shared objects (macOS/Linux)
- `*.dll` - Dynamic libraries (Windows)
- `symbols.rds` - Symbol cache

These files are automatically generated during compilation and will be different for:
- Intel Macs (x86_64)
- Apple Silicon Macs (arm64)
- Linux systems
- Windows systems

### .gitignore Recommendations

Make sure your `.gitignore` includes:
```
# Compiled files
src/*.o
src/*.so
src/*.dll
src/symbols.rds
src/*.rds

# R build artifacts
*.tar.gz
*.Rcheck/
.Rproj.user/
```

---

## System Information

**Architecture:** ARM64 (Apple Silicon)
**R Version:** 4.5.2
**Build Date:** December 1, 2025
**Compiler:** Apple clang version 17.0.0

---

## Troubleshooting

### If you get "architecture mismatch" error again:

1. Clean compiled files:
   ```bash
   cd src/
   rm -f *.o *.so symbols.rds
   ```

2. Rebuild:
   ```bash
   cd ..
   R CMD INSTALL --preclean .
   ```

### If compilation fails:

1. Check Rcpp and RcppArmadillo are installed:
   ```r
   install.packages(c("Rcpp", "RcppArmadillo"))
   ```

2. Check you have Xcode Command Line Tools:
   ```bash
   xcode-select --install
   ```

3. Verify R is configured for your architecture:
   ```bash
   R CMD config CC
   R CMD config CXX
   ```

---

## Package Testing

Now that the package is compiled correctly, you can run the test suites:

### Quick test (2-3 minutes):
```r
source("quick_test.R")
```

### Comprehensive test (10-15 minutes):
```r
source("comprehensive_tutorial.R")
```

Both test files will now work correctly on your ARM64 Mac.

---

**Status:** ✓ Ready for testing and development
