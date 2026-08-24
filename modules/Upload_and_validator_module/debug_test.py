"""Throwaway diagnostic script — walks through each validation step
individually, printing progress and catching everything (including
BaseException, which a normal 'except Exception' would miss — e.g. an
uncaught SystemExit would otherwise terminate the process with zero
output and no traceback, which matches what we're seeing)."""

import sys
import time

print("STEP 0: starting", flush=True)

try:
    from metabolomics_validator import InputSet, MetabolomicsValidator
    print("STEP 1: import ok", flush=True)

    inputs = InputSet(
        mzml_files=["sample_data/valid_files/Kgali_MeOH_pos_MCL_R1_1a.mzML"],
        metadata_file="sample_data/valid_files/Practice_data_MeOH_POS_metadata.csv",
        data_matrix_file="sample_data/valid_files/DataMatrix_Example.csv",
    )
    print("STEP 2: InputSet built", flush=True)

    validator = MetabolomicsValidator(inputs)
    print("STEP 3: validator constructed", flush=True)

    t0 = time.time()
    print("STEP 4: calling check_upload()...", flush=True)
    upload_ok = validator.check_upload()
    print(f"STEP 5: check_upload() done in {time.time()-t0:.2f}s -> {upload_ok}", flush=True)

    t0 = time.time()
    print("STEP 6: calling check_file_types()...", flush=True)
    types_ok = validator.check_file_types()
    print(f"STEP 7: check_file_types() done in {time.time()-t0:.2f}s -> {types_ok}", flush=True)

    t0 = time.time()
    print("STEP 8: calling _check_mzml_files() (the 58MB file — may take a while)...", flush=True)
    mzml_ok = validator._check_mzml_files()
    print(f"STEP 9: _check_mzml_files() done in {time.time()-t0:.2f}s -> {mzml_ok}", flush=True)

    t0 = time.time()
    print("STEP 10: calling _check_metadata()...", flush=True)
    meta_ok = validator._check_metadata()
    print(f"STEP 11: _check_metadata() done in {time.time()-t0:.2f}s -> {meta_ok}", flush=True)

    t0 = time.time()
    print("STEP 12: calling _check_data_matrix()...", flush=True)
    matrix_ok = validator._check_data_matrix()
    print(f"STEP 13: _check_data_matrix() done in {time.time()-t0:.2f}s -> {matrix_ok}", flush=True)

    print("STEP 14: all integrity checks completed successfully", flush=True)

    t0 = time.time()
    print("STEP 14a: calling check_compatibility()...", flush=True)
    compat_ok = validator.check_compatibility()
    print(f"STEP 14b: check_compatibility() done in {time.time()-t0:.2f}s -> {compat_ok}", flush=True)

    for issue in validator.issues:
        print(f"  ISSUE: [{issue.severity.value}] {issue.step.value} | {issue.file}: {issue.message}", flush=True)

    print("STEP 15: now calling validator.run() fresh, exactly like main() does...", flush=True)
    t0 = time.time()
    report = validator.run()
    print(f"STEP 16: run() done in {time.time()-t0:.2f}s, passed={report.passed}", flush=True)

    print("STEP 17: building report.to_text()...", flush=True)
    text = report.to_text()
    print(f"STEP 18: to_text() built ({len(text)} chars), about to print it", flush=True)
    print(text, flush=True)
    print("STEP 19: printed the report successfully", flush=True)

except BaseException as e:
    print(f"CAUGHT: {type(e).__module__}.{type(e).__name__}: {e!r}", flush=True)
    import traceback
    traceback.print_exc()
    if isinstance(e, SystemExit):
        print(f"SystemExit code was: {e.code!r}", flush=True)

print("STEP 20: script reached the very end", flush=True)