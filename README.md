bio-info-gpu
============

DNA/RNA to Amino Acids Translation using GPU CUDA

The objective of this program is to pratice some CUDA skills and write a simple a simple program using BioInformatics concepts.

The program receive an DNA or RNA Sequence and write the correspondent Amino Acid chain on output.

It breaks the input sequence in chunks and translate it in the GPU, brings back to CPU and write in the output file.

In order to run, you must import the code in NSight Visual Studio Edition and setup the NVidia SDK 5 (not tested in higher version, but it may work).

After compiled and builded, the command lines are:
program_name [path_to_sequence_file] [opt output_type] [opt path_to_output_file] [opt debug]

Where ouput_type assumes the following values: "compact" or "complete"

In complete output type each amino acid is shown on output with three letters (e.g: 'MET' for Methionine).
While in the compact output type shown just one letter (e.g: 'M' for Methionine).

Example:
Sequence: AAACCCTTTGGG
Using compact output type: KPFG
Using complete output type: LysProPheGly
