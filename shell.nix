{ pkgs ? import <nixpkgs> {}}:

with pkgs; mkShell {
  buildInputs = [
    clang-tools
  ] ++ lib.optionals stdenv.isDarwin [
    darwin.apple_sdk.frameworks.CoreFoundation
    darwin.apple_sdk.frameworks.Security
    darwin.apple_sdk.frameworks.SystemConfiguration
  ];
}
