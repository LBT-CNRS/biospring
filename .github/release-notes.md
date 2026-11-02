## Prebuilt executables

`biospring` and its tools, built with MDDriver, FreeSASA and OpenMP. Every non-system library is bundled, so nothing else needs to be installed. Unpack, then run `bin/biospring`.

| Archive | Platform | Minimum system |
|---|---|---|
| `*-macos-arm64.tar.gz` | Apple Silicon | macOS 14 |
| `*-macos-x86_64.tar.gz` | Intel Mac | macOS 15 |
| `*-linux-x86_64.tar.gz` | Linux, x86_64 | glibc 2.35 (Ubuntu 22.04, Debian 12 or newer) |
| `*-linux-aarch64.tar.gz` | Linux, ARM64 | glibc 2.35 (Ubuntu 22.04, Debian 12 or newer) |

Verify a download with the matching `.sha256` file (`shasum -a 256 -c` on macOS, `sha256sum -c` on Linux).

**macOS:** the binaries carry an ad-hoc signature and are not notarized, so a browser download is blocked by Gatekeeper, once per file. After unpacking, clear the quarantine flag for the whole folder in one go:

```
xattr -dr com.apple.quarantine biospring-*-macos-*/
```

Downloading with `curl -L -O <url>` does not set the flag at all.

**Windows:** there is no native build. The Linux x86_64 archive is expected to run under WSL2 (not tested).

---

