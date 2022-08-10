# Add inputs and outputs from these tool invocations to the build variables
CPP_SRCS += \
./src/fasta2Fmi.cpp

OBJS += \
./src/fasta2Fmi.o

CPP_DEPS += \
./src/fasta2Fmi.d


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ./src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_BZIP2=0 -DSEQAN_HAS_ZLIB=1 -O3 -pedantic -Wall -c -fmessage-length=0 -W -Wno-long-long -Wno-variadic-macros -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
