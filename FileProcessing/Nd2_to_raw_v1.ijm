/*
 * Batch-conversion from .nd2 to raw using the Raw-Yaml Exporter from jungmannlab (https://github.com/jungmannlab/imagej-raw-yaml-export)
 * Relatively fast due to virtual loading of .nd2 file in batch-mode
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".nd2") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	print("Processing: " + input + File.separator + file);
	//setBatchMode(true);
	run("Bio-Formats", "open=[" + input + "/" + file +"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
	//title_long = getTitle();
	title = split(getTitle(), "/");
	print(title[1]);
	rename(title[1]);
	run("raw-yaml Exporter", "save=[" + output + "/" + title[1] + ".raw"+"]");
	print("Saving to: " + output);
	close("*");
	}
