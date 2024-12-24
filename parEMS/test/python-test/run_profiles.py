import subprocess, random


def generate_test_file(file_name, num_sequences, length_range, seed=2333):
    random.seed(seed)
    low, high = length_range
    random_range = lambda low, high: lambda: random.randint(low, high)
    with open(file_name, 'w') as file:
        for _ in range(num_sequences):
            sequence = ''.join(random.choice('ATGC') for _ in range(random_range(low, high)()))
            file.write(sequence + '\n')
    return file_name


def run_cli(binary_path, string):
    command = f"{binary_path} {string}".split()
    result = subprocess.run(command, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


def extract_usage(stdout):
    outputs = {}
    for line in stdout.splitlines():
        line = line.strip()
        if line.startswith("l = "):
            # read config that are
            for part in line.split(", "):
                kind, value = part.split(" = ")
                outputs[kind] = int(value)
        if "TOTAL" in line:
            if "Time" in line:
                outputs["seconds"] = float(line.split(" ")[-2])
            else:
                outputs["gigabytes"] = float(line.split(" ")[-2])
    return outputs


def test_configurations():
    file_path = generate_test_file('@CLI_BINARY_PATH@/testfile', 100, (20, 35), 2333)

    method_configs = [
        # "1",  # EMS1
        # "2",  # EMS2
        # "4",  # RecursEMS
        "5",  # ParEMS
    ]

    ld_configs = [
        (l, d) for l in range(2, 8) for d in range(l)
    ]

    test_configs = [
        f"-s {m} -l {l} -d {d} -t {t} {file_path}"
        for t in range(10) for l, d in ld_configs for m in method_configs
    ]

    # test_configs = [
    #     f"-s 5 -l 4 -d 0 -t 1 {file_path}"
    # ]

    with open("test_results.csv", "w") as f:
        f.write("l, d, t, secs, GB\n")
        for config in test_configs:
            returncode, stdout, stderr = run_cli("@CLI_BINARY_PATH@/parEMS", config)
            if returncode != 0:
                print(f"Test failed with config {config}: {stderr}")
                continue
            # assert "expected_output" in stdout, f"Output mismatch with config {config}"
            usage = extract_usage(stdout)
            # write usage to file
            f.write(f"{usage['l']}, {usage['d']}, "
                    f"{usage['t']}, {usage['seconds']:0.5f}, "
                    f"{usage['gigabytes']:0.5f}\n")


if __name__ == "__main__":
    test_configurations()
