# Wrapper Makefile for CMake
.PHONY: all build clean rebuild run
BUILD_DIR ?= build
all: build
configure:
	@mkdir -p $(BUILD_DIR)
	cmake -S . -B $(BUILD_DIR)
build: configure
	cmake --build $(BUILD_DIR) -- -j
clean:
	@rm -rf $(BUILD_DIR)
rebuild: clean build
run-plots: build
	@./$(BUILD_DIR)/ddis_plots_q2 input.root
run-skim: build
	@./$(BUILD_DIR)/ddis_skim_q2 data/filelist.txt skim_q2.root 100
