FC ?= gfortran
FFLAGS ?= -std=legacy -ffixed-form -ffixed-line-length-none -O2

.PHONY: all run clean

all: eaf_general_balance

eaf_general_balance: eaf_general_balance.f
	$(FC) $(FFLAGS) $< -o $@

run: eaf_general_balance
	./eaf_general_balance

clean:
	rm -f eaf_general_balance eaf_general_balance.exe
	rm -f general_eaf_output.txt
