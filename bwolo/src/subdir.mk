# Add inputs and outputs from these tool invocations to the build variables
CPP_SRCS += \
./src/AlignmentInFmi.cpp \
./src/CandidateLog.cpp \
./src/ExtendRight.cpp \
./src/PatternFiltrer.cpp \
./src/PatternSegment.cpp \
./src/SegmentErrBrowse.cpp \
./src/interfaceMainFonction.cpp \
./src/bwolo.cpp \
./src/seqAlign.cpp

OBJS += \
./src/AlignmentInFmi.o \
./src/CandidateLog.o \
./src/ExtendRight.o \
./src/PatternFiltrer.o \
./src/PatternSegment.o \
./src/SegmentErrBrowse.o \
./src/interfaceMainFonction.o \
./src/bwolo.o \
./src/seqAlign.o

CPP_DEPS += \
./src/AlignmentInFmi.d \
./src/CandidateLog.d \
./src/ExtendRight.d \
./src/PatternFiltrer.d \
./src/PatternSegment.d \
./src/SegmentErrBrowse.d \
./src/interfaceMainFonction.d \
./src/bwolo.d \
./src/seqAlign.d

STD=c++11

# Each subdirectory must supply rules for building sources it contributes
src/%.o: ./src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -isystem ../seqan/include -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_BZIP2=0 -DSEQAN_HAS_ZLIB=1 -O3 -pedantic -Wall -c -fmessage-length=0 -W -Wno-long-long -Wno-variadic-macros -std=$(STD) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
