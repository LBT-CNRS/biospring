class Biospring < Formula
  desc "Coarse-grained interactive molecular simulation (with MDDriver, FreeSASA, OpenMP)"
  homepage "https://github.com/LBT-CNRS/biospring"
  version "1.1.2"
  license "CECILL-C"

  # Prebuilt, relocatable archives produced by .github/workflows/release.yml.
  # Update version and all four sha256 values (the release ships *.sha256 files) for each release.
  on_macos do
    on_arm do
      url "https://github.com/LBT-CNRS/biospring/releases/download/v#{version}/biospring-v#{version}-macos-arm64.tar.gz"
      sha256 "0000000000000000000000000000000000000000000000000000000000000000"
    end
    on_intel do
      url "https://github.com/LBT-CNRS/biospring/releases/download/v#{version}/biospring-v#{version}-macos-x86_64.tar.gz"
      sha256 "0000000000000000000000000000000000000000000000000000000000000000"
    end
  end
  on_linux do
    on_arm do
      url "https://github.com/LBT-CNRS/biospring/releases/download/v#{version}/biospring-v#{version}-linux-aarch64.tar.gz"
      sha256 "0000000000000000000000000000000000000000000000000000000000000000"
    end
    on_intel do
      url "https://github.com/LBT-CNRS/biospring/releases/download/v#{version}/biospring-v#{version}-linux-x86_64.tar.gz"
      sha256 "0000000000000000000000000000000000000000000000000000000000000000"
    end
  end

  def install
    prefix.install Dir["*"]
  end

  test do
    assert_match "biospring", shell_output("#{bin}/biospring --version 2>&1")
  end
end
