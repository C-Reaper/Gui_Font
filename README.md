# Project README

## Overview
This project is a simple application that renders "Hello World!" using different font styles in various platforms (Linux, Windows, Wine, and Web). It demonstrates the use of custom fonts and rendering functions in a cross-platform manner.

## Features
- Renders text "Hello World!" in three different font styles.
- Supports Linux, Windows, Wine, and Web platforms.
- Uses custom font loading and rendering functions.

## Project Structure
- `build/`: Stores the executable files produced by `Main.c`.
- `src/`: Contains the source code.
  - `Main.c`: Entry point of the application.
  - `*.h`: Standalone header-based C-files without corresponding `.c` files that implement them.
- `Makefile.linux`: Linux build configuration.
- `Makefile.windows`: Windows build configuration.
- `Makefile.wine`: Wine build configuration for cross-compiling to Windows on Linux.
- `Makefile.web`: Emscripten build configuration for WebAssembly.

## Prerequisites
- C/C++ Compiler and Debugger (GCC, Clang)
- Make utility
- Standard development tools
- Libraries needed in specific projects:
  - **Linux**: X11, PNG, JPEG
  - **Windows**: WINAPI (for Windows API), GDI32, USER32, WINMM
  - **Wine**: No additional libraries required since it uses the same APIs as Wine.
  - **Web**: Emscripten SDK

## Build & Run
### Linux
To build and run on Linux:
```sh
cd <Project>
make -f Makefile.linux all
make -f Makefile.linux exe
```

### Windows
To build and run on Windows:
```sh
cd <Project>
make -f Makefile.windows all
./build/Main.exe
```

### Wine
To cross-compile for Windows using Wine on Linux:
```sh
cd <Project>
make -f Makefile.wine all
wine ./build/Main.exe
```

### Web
To build and run in the browser (WebAssembly):
```sh
cd <Project>
make -f Makefile.web all
./build/index.html
```

This setup ensures that the project can be built and executed across different platforms using standard tools and libraries.