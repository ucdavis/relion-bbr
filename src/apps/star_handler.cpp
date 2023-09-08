/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <src/image.h>
#include <src/metadata_table.h>
#include <src/filename.h>
#include <src/time.h>
#include <src/jaz/obs_model.h>
#include <src/pipeline_jobs.h>
#include <src/symmetries.h>
#include <cmath>

class star_handler_parameters
{
	public:

	FileName fn_in, tablename_in, fn_out, fn_compare, tablename_compare,
	         fn_label1, fn_label2, fn_label3,
	         select_label, select_str_label, discard_label,
	         fn_check, fn_operate, fn_operate2, fn_operate3, fn_set;

	std::string remove_col_label, add_col_label, add_col_value, add_col_from, hist_col_label, select_include_str, select_exclude_str;
	RFLOAT eps, select_minval, select_maxval, multiply_by, add_to, center_X, center_Y, center_Z, hist_min, hist_max;
	bool do_ignore_optics, do_combine, do_split, do_center, do_random_order, show_frac, show_cumulative, do_discard,do_RoT,do_BBR,do_remove_nan;
	long int nr_split, size_split, nr_bin, random_seed;
	RFLOAT discard_sigma, duplicate_threshold, extract_angpix, cl_angpix,rOt,tIlt,pSi;
	RFLOAT BBR_X,BBR_Y,BBR_Z,micrograph_Xsize,micrograph_Ysize;
	bool handerness,handerness_negative,do_invert_handerness,do_microgrpah_drop,ROUND_SHIFT;
	bool do_write_bild;
	RFLOAT bild_angle_step,bild_radius;
	FileName fn_sym;
	ObservationModel obsModel;
	// I/O Parser
	IOParser parser;
	Matrix1D<RFLOAT> Posi;
	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
		fn_in = parser.getOption("--i", "Input STAR file");
		fn_out = parser.getOption("--o", "Output STAR file", "out.star");
		do_ignore_optics = parser.checkOption("--ignore_optics", "Provide this option for relion-3.0 functionality, without optics groups");
		cl_angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstrom, for when ignoring the optics groups in the input star file", "1."));
		tablename_in = parser.getOption("--i_tablename", "If ignoring optics, then read table with this name", "");

		int compare_section = parser.addSection("Compare options");
		fn_compare = parser.getOption("--compare", "STAR file name to compare the input STAR file with", "");
		fn_label1 = parser.getOption("--label1", "1st metadata label for the comparison (may be string, int or RFLOAT)", "");
		fn_label2 = parser.getOption("--label2", "2nd metadata label for the comparison (RFLOAT only) for 2D/3D-distance)", "");
		fn_label3 = parser.getOption("--label3", "3rd metadata label for the comparison (RFLOAT only) for 3D-distance)", "");
		eps = textToFloat(parser.getOption("--max_dist", "Maximum distance to consider a match (for int and RFLOAT only)", "0."));
		
		////
		int Rot_star_with_rtp_section= parser.addSection("Rotate star to given alignment");
		do_RoT = parser.checkOption("--do_rotation", "Provide this option to apply rOt, tIlt, pSi to the starfile.");
		rOt = textToFloat(parser.getOption("--rOt", "rot angle", "0."));
		tIlt = textToFloat(parser.getOption("--tIlt", "rot angle", "0."));
		pSi = textToFloat(parser.getOption("--pSi", "rot angle", "0."));
		////
		int BBR_section= parser.addSection("!!!!Block-based reconstruction");
		do_BBR = parser.checkOption("--do_BBR", "Provide this option to produce star files for block-based reconstruction.");
		fn_sym = parser.getOption("--sym", "Symmetry of input STAR file","");
		BBR_X = textToFloat(parser.getOption("--BBR_X", "X position in pixel", "0."));
		BBR_Y = textToFloat(parser.getOption("--BBR_Y", "Y position in pixel", "0."));
		BBR_Z = textToFloat(parser.getOption("--BBR_Z", "Z position in pixel", "0."));
		handerness = parser.checkOption("--correct_handerness", "Suppose your handerness was incorrect.");
		handerness_negative = parser.checkOption("--handerness_negative", "Suppose your handerness was inverted. If both options are missing, suppose no change to defocus.");
		micrograph_Xsize = textToFloat(parser.getOption("--micrograph_xsize", "xsize of the micrographs", "-1."));
		micrograph_Ysize = textToFloat(parser.getOption("--micrograph_ysize", "xsize of the micrographs", "-1."));
		ROUND_SHIFT = parser.checkOption("--round_shift", "Output only the Interger part of the X/Y shift of particles.");
		int SP_BBR_section= parser.addSection("!!!!Special option of Block-based reconstruction");
		do_invert_handerness = parser.checkOption("--do_invert_BBR_handerness", "Suppose your handerness was correct.");
		
		int bild_section=parser.addSection("Write bild options");
		do_write_bild=parser.checkOption("--write_bild", "Generate bild from input star file.");
		bild_angle_step=textToFloat(parser.getOption("--bild_angle_step", "Angle step of bild", "1.875"));
		bild_radius=textToFloat(parser.getOption("--bild_radius", "Radius of inner sphere of bild", "128.0"));
		
		int subset_section = parser.addSection("Select options");
		select_label = parser.getOption("--select", "Metadata label (number) to base output selection on (e.g. rlnCtfFigureOfMerit)", "");
		select_minval = textToFloat(parser.getOption("--minval", "Minimum acceptable value for this label (inclusive)", "-99999999."));
		select_maxval = textToFloat(parser.getOption("--maxval", "Maximum acceptable value for this label (inclusive)", "99999999."));
		select_str_label = parser.getOption("--select_by_str", "Metadata label (string) to base output selection on (e.g. rlnMicrographname)", "");
		select_include_str = parser.getOption("--select_include", "select rows that contains this string in --select_by_str ", "");
		select_exclude_str = parser.getOption("--select_exclude", "exclude rows that contains this string in --select_by_str ", "");

		int discard_section = parser.addSection("Discard based on image statistics options");
		do_discard = parser.checkOption("--discard_on_stats", "Discard images if their average/stddev deviates too many sigma from the ensemble average");
		do_remove_nan = parser.checkOption("--remove_nan", "remove images that has nan in them from star file. Should combine with --discard_on_stats");
		discard_label = parser.getOption("--discard_label", "MetaDataLabel that points to the images to be used for discarding based on statistics", "rlnImageName");
		discard_sigma = textToFloat(parser.getOption("--discard_sigma", "Discard images with average or stddev values that lie this many sigma away from the ensemble average", "4."));

		int combine_section = parser.addSection("Combine options");
		do_combine = parser.checkOption("--combine", "Combine input STAR files (multiple individual filenames, all within double-quotes after --i)");
		fn_check = parser.getOption("--check_duplicates", "MetaDataLabel (for a string only!) to check for duplicates, e.g. rlnImageName", "");

		int split_section = parser.addSection("Split options");
		do_split = parser.checkOption("--split", "Split the input STAR file into one or more smaller output STAR files");
		do_random_order = parser.checkOption("--random_order", "Perform splits on randomised order of the input STAR file");
		random_seed = textToInteger(parser.getOption("--random_seed", "Random seed for randomisation.", "-1"));
		nr_split = textToInteger(parser.getOption("--nr_split", "Split into this many equal-sized STAR files", "-1"));
		size_split = textToLongLong(parser.getOption("--size_split", "AND/OR split into subsets of this many lines", "-1"));

		int operate_section = parser.addSection("Operate options");
		fn_operate = parser.getOption("--operate", "Operate on this metadata label", "");
		fn_operate2 = parser.getOption("--operate2", "Operate also on this metadata label", "");
		fn_operate3 = parser.getOption("--operate3", "Operate also on this metadata label", "");
		fn_set = parser.getOption("--set_to", "Set all the values for the --operate label(s) to this value", "");
		multiply_by = textToFloat(parser.getOption("--multiply_by", "Multiply all the values for the --operate label(s) by this value", "1."));
		add_to = textToFloat(parser.getOption("--add_to", "Add this value to all the values for the --operate label(s)", "0."));

		int center_section = parser.addSection("Center options");
		do_center = parser.checkOption("--center", "Perform centering of particles according to a position in the reference.");
		center_X = textToFloat(parser.getOption("--center_X", "X-coordinate in the reference to center particles on (in pix)", "0."));
		center_Y = textToFloat(parser.getOption("--center_Y", "Y-coordinate in the reference to center particles on (in pix)", "0."));
		center_Z = textToFloat(parser.getOption("--center_Z", "Z-coordinate in the reference to center particles on (in pix)", "0."));

		int column_section = parser.addSection("Column options");
		remove_col_label = parser.getOption("--remove_column", "Remove the column with this metadata label from the input STAR file.", "");
		add_col_label = parser.getOption("--add_column", "Add a column with this metadata label from the input STAR file.", "");
		add_col_value = parser.getOption("--add_column_value", "Set this value in all rows for the added column", "");
		add_col_from = parser.getOption("--copy_column_from", "Copy values in this column to the added column", "");
		hist_col_label = parser.getOption("--hist_column", "Calculate histogram of values in the column with this metadata label", "");
		show_frac = parser.checkOption("--in_percent", "Show a histogram in percent (need --hist_column)");
		show_cumulative = parser.checkOption("--show_cumulative", "Show a histogram of cumulative distribution (need --hist_column)");
		nr_bin = textToInteger(parser.getOption("--hist_bins", "Number of bins for the histogram. By default, determined automatically by Freedman–Diaconis rule.", "-1"));
		hist_min = textToFloat(parser.getOption("--hist_min", "Minimum value for the histogram (needs --hist_bins)", "-inf"));
		hist_max = textToFloat(parser.getOption("--hist_max", "Maximum value for the histogram (needs --hist_bins)", "inf"));

		int duplicate_section = parser.addSection("Duplicate removal");
		duplicate_threshold = textToFloat(parser.getOption("--remove_duplicates","Remove duplicated particles within this distance [Angstrom]. Negative values disable this.", "-1"));
		extract_angpix = textToFloat(parser.getOption("--image_angpix", "For down-sampled particles, specify the pixel size [A/pix] of the original images used in the Extract job", "-1"));

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line, exiting...");
	}

	void run()
	{
		int c = 0;
		if (fn_compare != "") c++;
		if (select_label != "") c++;
		if (select_str_label != "") c++;
		if (do_discard) c++;
		if (do_combine) c++;
		if (do_split) c++;
		if (fn_operate != "") c++;
		if (do_center) c++;
		if (remove_col_label != "") c++;
		if (add_col_label != "") c++;
		if (hist_col_label != "") c++;
		if (duplicate_threshold > 0) c++;
		if (c != 1)
		{
			MetaDataTable MD;
			read_check_ignore_optics(MD, fn_in);
			write_check_ignore_optics(MD, fn_out, MD.getName());
			//REPORT_ERROR("ERROR: specify (only and at least) one of the following options: --compare, --select, --select_by_str, --combine, --split, --operate, --center, --remove_column, --add_column, --hist_column or --remove_duplicates.");
		}

		if (fn_out == "" && hist_col_label == "")
			REPORT_ERROR("ERROR: specify the output file name (--o)");

		if (fn_compare != "") compare();
		if (select_label != "") select();
		if (select_str_label != "") select_by_str();
		if (do_discard) discard_on_image_stats();
		////
		if(do_write_bild) do_bild();
		////
		if(do_RoT) do_rot_star();
		////
		////
		if(do_BBR){
			if(handerness && handerness_negative){
				REPORT_ERROR("ERROR: --correct_handerness and --handerness_negative cannot be used at the same time.");
			}
			Posi.resize(3);
			Posi(0)=BBR_X;
			Posi(1)=BBR_Y;
			Posi(2)=BBR_Z;
			do_microgrpah_drop=false;
			if(micrograph_Xsize>0 && micrograph_Ysize>0){
				do_microgrpah_drop=true;
			}
			do_block_based();
		} 
		////
		if(do_invert_handerness) invert_BBR_handerness();
		if (do_combine) combine();
		if (do_split) split();
		if (fn_operate != "") operate();
		if (do_center) center();
		if (remove_col_label!= "") remove_column();
		if (add_col_label!= "") add_column();
		if (hist_col_label != "") hist_column();
		if (duplicate_threshold > 0) remove_duplicate();

		std::cout << " Done!" << std::endl;
	}

	void read_check_ignore_optics(MetaDataTable &MD, FileName fn, std::string tablename = "discover")
	{
		if (do_ignore_optics)
		{
			MD.read(fn, tablename);
		}
		else
		{
			ObservationModel::loadSafely(fn, obsModel, MD, tablename, 1, false);
			if (obsModel.opticsMdt.numberOfObjects() == 0)
			{
				std::cerr << " + WARNGING: could not read optics groups table, proceeding without it ..." << std::endl;
				MD.read(fn, tablename);
				do_ignore_optics = true;
			}
		}
	}

	void write_check_ignore_optics(MetaDataTable &MD, FileName fn, std::string tablename)
	{
		if (do_ignore_optics) MD.write(fn);
		else obsModel.save(MD, fn, tablename);
	}

	void compare()
	{
	   	MetaDataTable MD1, MD2, MDonly1, MDonly2, MDboth;
		EMDLabel label1, label2, label3;

		// Read in the observationModel
		read_check_ignore_optics(MD2, fn_compare);
		// read_check_ignore_optics() overwrites the member variable obsModel (BAD DESIGN!)
		// so we have to back up.
		ObservationModel obsModelCompare = obsModel;
		read_check_ignore_optics(MD1, fn_in);

		label1 = EMDL::str2Label(fn_label1);
		label2 = (fn_label2 == "") ? EMDL_UNDEFINED : EMDL::str2Label(fn_label2);
		label3 = (fn_label3 == "") ? EMDL_UNDEFINED : EMDL::str2Label(fn_label3);

		compareMetaDataTable(MD1, MD2, MDboth, MDonly1, MDonly2, label1, eps, label2, label3);

		std::cout << MDboth.numberOfObjects()  << " entries occur in both input STAR files." << std::endl;
		std::cout << MDonly1.numberOfObjects() << " entries occur only in the 1st input STAR file." << std::endl;
		std::cout << MDonly2.numberOfObjects() << " entries occur only in the 2nd input STAR file." << std::endl;

		write_check_ignore_optics(MDboth, fn_out.insertBeforeExtension("_both"), MD1.getName());
		std::cout << " Written: " << fn_out.insertBeforeExtension("_both") << std::endl;
		write_check_ignore_optics(MDonly1, fn_out.insertBeforeExtension("_only1"), MD1.getName());
		std::cout << " Written: " << fn_out.insertBeforeExtension("_only1") << std::endl;
		// Use MD2's optics group for MDonly2.
		obsModel = obsModelCompare;
		write_check_ignore_optics(MDonly2, fn_out.insertBeforeExtension("_only2"), MD1.getName());
		std::cout << " Written: " << fn_out.insertBeforeExtension("_only2") << std::endl;
	}

	void select()
	{
		MetaDataTable MDin, MDout;

		read_check_ignore_optics(MDin, fn_in);

		MDout = subsetMetaDataTable(MDin, EMDL::str2Label(select_label), select_minval, select_maxval);

		write_check_ignore_optics(MDout, fn_out, MDin.getName());
		std::cout << " Written: " << fn_out << " with " << MDout.numberOfObjects() << " item(s)" << std::endl;
	}

	void select_by_str()
	{
		int c = 0;
		if (select_include_str != "") c++;
		if (select_exclude_str != "") c++;
		if (c != 1)
			REPORT_ERROR("You must specify only and at least one of --select_include and --select_exclude");

		MetaDataTable MDin, MDout;

		read_check_ignore_optics(MDin, fn_in);

		if (select_include_str != "")
			MDout = subsetMetaDataTable(MDin, EMDL::str2Label(select_str_label), select_include_str, false);
		else
			MDout = subsetMetaDataTable(MDin, EMDL::str2Label(select_str_label), select_exclude_str, true);

		write_check_ignore_optics(MDout, fn_out, MDin.getName());
		std::cout << " Written: " << fn_out << std::endl;

	}

	void discard_on_image_stats()
	{
		MetaDataTable MDin, MDout;
		read_check_ignore_optics(MDin, fn_in);

		std::cout << " Calculating average and stddev for all images ... " << std::endl;
		time_config();
		init_progress_bar(MDin.numberOfObjects());


   		RFLOAT sum_avg = 0.;
		RFLOAT sum2_avg = 0.;
		RFLOAT sum_stddev = 0.;
		RFLOAT sum2_stddev = 0.;
		RFLOAT sum_n = 0.;
		RFLOAT hasnan=0.;
		std::vector<RFLOAT> avgs, stddevs,nan_discard;
		long int ii = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			Image<RFLOAT> img;
			FileName fn_img;
			RFLOAT avg, stddev, minval, maxval;
			MDin.getValue(EMDL::str2Label(discard_label), fn_img);
			img.read(fn_img);
			hasnan=0.0;
			if(do_remove_nan){
				img().setXmippOrigin();
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
				{
					if (std::isnan(DIRECT_A3D_ELEM(img(), k, i, j)) || std::isinf(DIRECT_A3D_ELEM(img(), k, i, j))){
						DIRECT_A3D_ELEM(img(), k, i, j) = 0.0;
						hasnan+=1.0;
					}
				}
			}
			img().computeStats(avg, stddev, minval, maxval);
			sum_avg += avg;
			sum2_avg += avg * avg;
			sum_stddev += stddev;
			sum2_stddev += stddev * stddev;
			sum_n += 1.;
			avgs.push_back(avg);
			stddevs.push_back(stddev);
			nan_discard.push_back(hasnan);
			ii++;
			if (ii%100 == 0) progress_bar(ii);
		}

		progress_bar(MDin.numberOfObjects());

		sum_avg /= sum_n;
		sum_stddev /= sum_n;
		sum2_avg = sqrt(sum2_avg/sum_n - sum_avg*sum_avg);
		sum2_stddev = sqrt(sum2_stddev/sum_n - sum_stddev*sum_stddev);

		std::cout << " [ Average , stddev ] of the average Image value = [ " << sum_avg<< " , " << sum2_avg << " ] " << std::endl;
		std::cout << " [ Average , stddev ] of the stddev  Image value = [ " << sum_stddev<< " , " << sum2_stddev << " ] "  << std::endl;

		long int i = 0, nr_discard = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			if (avgs[i] > sum_avg - discard_sigma * sum2_avg &&
			    avgs[i] < sum_avg + discard_sigma * sum2_avg &&
			    stddevs[i] > sum_stddev - discard_sigma * sum2_stddev &&
			    stddevs[i] < sum_stddev + discard_sigma * sum2_stddev && nan_discard[i]<0.5)
			{
				MDout.addObject(MDin.getObject(current_object));
			}
			else
			{
				nr_discard++;
				if(nan_discard[i]>0.5){
					FileName fn_img;
					MDin.getValue(EMDL_IMAGE_NAME, fn_img);
					std::cout<<"Has nan on "<<fn_img<<std::endl;
				}
			}
			i++;
		}

		std::cout << " Discarded " << nr_discard << " Images because of too large or too small average/stddev values " << std::endl;

		write_check_ignore_optics(MDout, fn_out, MDin.getName());
		std::cout << " Written: " << fn_out << std::endl;
	}

	void combine()
	{

		std::vector<FileName> fns_in;
		std::vector<std::string> words;
		tokenize(fn_in, words);
		for (int iword = 0; iword < words.size(); iword++)
		{
			FileName fnt = words[iword];
			fnt.globFiles(fns_in, false);
		}

		MetaDataTable MDin, MDout;
		std::vector<MetaDataTable> MDsin, MDoptics;
		std::vector<ObservationModel> obsModels;
		// Read the first table into the global obsModel
		read_check_ignore_optics(MDin, fns_in[0]);
		MDsin.push_back(MDin);
		// Read all the rest of the tables into local obsModels
		for (int i = 1; i < fns_in.size(); i++)
		{
			ObservationModel myobsModel;
			if (do_ignore_optics) MDin.read(fns_in[i], tablename_in);
			else ObservationModel::loadSafely(fns_in[i], myobsModel, MDin, "discover", 1);
			MDsin.push_back(MDin);
			obsModels.push_back(myobsModel);
		}

		// Combine optics groups with the same EMDL_IMAGE_OPTICS_GROUP_NAME, make new ones for those with a different name
		if (!do_ignore_optics)
		{
			std::vector<std::string> optics_group_uniq_names;

			// Initialise optics_group_uniq_names with the first table
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
			{
				std::string myname;
				obsModel.opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP_NAME, myname);
				optics_group_uniq_names.push_back(myname);
			}

			// Now check uniqueness of the other tables
			for (int MDs_id = 1; MDs_id < fns_in.size(); MDs_id++)
			{
				const int obs_id = MDs_id - 1;

				std::vector<int> new_optics_groups;
				FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsin[MDs_id])
				{
					int tmp;
					MDsin[MDs_id].getValue(EMDL_IMAGE_OPTICS_GROUP, tmp);
					new_optics_groups.push_back(tmp);
				}

				MetaDataTable unique_opticsMdt;
				unique_opticsMdt.addMissingLabels(&obsModels[obs_id].opticsMdt);

				FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModels[obs_id].opticsMdt)
				{
					std::string myname;
					int my_optics_group;
					obsModels[obs_id].opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP_NAME, myname);
					obsModels[obs_id].opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, my_optics_group);

					// Check whether this name is unique
					bool is_uniq = true;
					int new_group;
					for (new_group = 0; new_group < optics_group_uniq_names.size(); new_group++)
					{
						if (optics_group_uniq_names[new_group] == myname)
						{
							is_uniq = false;
							break;
						}
					}
					new_group ++; // start counting of groups at 1, not 0!

					if (is_uniq)
					{
						std::cout << " + Adding new optics_group with name: " << myname << std::endl;

						optics_group_uniq_names.push_back(myname);
						// Add the line to the global obsModel
						obsModels[obs_id].opticsMdt.setValue(EMDL_IMAGE_OPTICS_GROUP, new_group);

						unique_opticsMdt.addObject();
						unique_opticsMdt.setObject(obsModels[obs_id].opticsMdt.getObject());
					}
					else
					{
						std::cout << " + Joining optics_groups with the same name: " << myname << std::endl;
						std::cerr << " + WARNING: if these are different data sets, you might want to rename optics groups instead of joining them!" << std::endl;
						std::cerr << " + WARNING: if so, manually edit the rlnOpticsGroupName column in the optics_groups table of your input STAR files." << std::endl;
					}

					if (my_optics_group != new_group)
					{
						std::cout << " + Renumbering group " << myname << " from " << my_optics_group << " to " << new_group << std::endl;
					}

					// Update the optics_group entry for all particles in the MDsin
					for (long int current_object2 = MDsin[MDs_id].firstObject();
					     current_object2 < MDsin[MDs_id].numberOfObjects() && current_object2 >= 0;
					     current_object2 = MDsin[MDs_id].nextObject())
					{
						int old_optics_group;
						MDsin[MDs_id].getValue(EMDL_IMAGE_OPTICS_GROUP, old_optics_group, current_object2);
						if (old_optics_group == my_optics_group)
							new_optics_groups[current_object2] = new_group;
					}
				}

				obsModels[obs_id].opticsMdt = unique_opticsMdt;

				FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsin[MDs_id])
				{
					MDsin[MDs_id].setValue(EMDL_IMAGE_OPTICS_GROUP, new_optics_groups[current_object]);

					// Also rename the rlnGroupName to not have groups overlapping from different optics groups
					std::string name;
					if (MDsin[MDs_id].getValue(EMDL_MLMODEL_GROUP_NAME, name))
					{
						name = "optics"+integerToString(new_optics_groups[current_object])+"_"+name;
						MDsin[MDs_id].setValue(EMDL_MLMODEL_GROUP_NAME, name);
					}
				}
			}

			// Make one vector for combination of the optics tables
			MDoptics.push_back(obsModel.opticsMdt);
			for (int i = 1; i < fns_in.size(); i++)
			{
				MDoptics.push_back(obsModels[i - 1].opticsMdt);
			}

			// Check if anisotropic magnification and/or beam_tilt are present in some optics groups, but not in others.
			// If so, initialise the others correctly
			bool has_beamtilt = false, has_not_beamtilt = false;
			bool has_anisomag = false, has_not_anisomag = false;
			bool has_odd_zernike = false, has_not_odd_zernike = false;
			bool has_even_zernike = false, has_not_even_zernike = false;
			bool has_ctf_premultiplied = false, has_not_ctf_premultiplied = false;
			for (int i = 0; i < fns_in.size(); i++)
			{
				if (MDoptics[i].containsLabel(EMDL_IMAGE_BEAMTILT_X) ||
					MDoptics[i].containsLabel(EMDL_IMAGE_BEAMTILT_Y))
				{
					has_beamtilt = true;
				}
				else
				{
					has_not_beamtilt = true;
				}
				if (MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_00) &&
				    MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_01) &&
				    MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_10) &&
				    MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_11))
				{
					has_anisomag = true;
				}
				else
				{
					has_not_anisomag = true;
				}
				if (MDoptics[i].containsLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS))
				{
					has_odd_zernike = true;
				}
				else
				{
					has_not_odd_zernike = true;
				}
				if (MDoptics[i].containsLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS))
				{
					has_even_zernike = true;
				}
				else
				{
					has_not_even_zernike = true;
				}
				if (MDoptics[i].containsLabel(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED))
				{
					has_ctf_premultiplied = true;
				}
				else
				{
					has_not_ctf_premultiplied = true;
				}
			}
#ifdef DEBUG
			printf("has_beamtilt = %d, has_not_beamtilt = %d, has_anisomag = %d, has_not_anisomag = %d, has_odd_zernike = %d, has_not_odd_zernike = %d, has_even_zernike = %d, has_not_even_zernike = %d, has_ctf_premultiplied = %d, has_not_ctf_premultiplied = %d\n", has_beamtilt, has_not_beamtilt, has_anisomag, has_not_anisomag, has_odd_zernike, has_not_odd_zernike, has_even_zernike, has_not_even_zernike, has_ctf_premultiplied, has_not_ctf_premultiplied);
#endif

			for (int i = 0; i < fns_in.size(); i++)
			{
				if (has_beamtilt && has_not_beamtilt)
				{
					if (!MDoptics[i].containsLabel(EMDL_IMAGE_BEAMTILT_X))
					{
						FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
						{
							MDoptics[i].setValue(EMDL_IMAGE_BEAMTILT_X, 0.);
						}
					}
					if (!MDoptics[i].containsLabel(EMDL_IMAGE_BEAMTILT_Y))
					{
						FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
						{
							MDoptics[i].setValue(EMDL_IMAGE_BEAMTILT_Y, 0.);
						}
					}
				}

				if (has_anisomag && has_not_anisomag)
				{
					if (!(MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_00) &&
					      MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_01) &&
					      MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_10) &&
					      MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_11)) )
					{
						FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
						{
							MDoptics[i].setValue(EMDL_IMAGE_MAG_MATRIX_00, 1.);
							MDoptics[i].setValue(EMDL_IMAGE_MAG_MATRIX_01, 0.);
							MDoptics[i].setValue(EMDL_IMAGE_MAG_MATRIX_10, 0.);
							MDoptics[i].setValue(EMDL_IMAGE_MAG_MATRIX_11, 1.);
						}
					}
				}

				if (has_odd_zernike && has_not_odd_zernike)
				{
					std::vector<RFLOAT> six_zeros(6, 0);
					if (!MDoptics[i].containsLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS))
					{
						FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
						{
							MDoptics[i].setValue(EMDL_IMAGE_ODD_ZERNIKE_COEFFS, six_zeros);
						}
					}
				}

				if (has_even_zernike && has_not_even_zernike)
				{
					std::vector<RFLOAT> nine_zeros(9, 0);
					if (!MDoptics[i].containsLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS))
					{
						FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
						{
							MDoptics[i].setValue(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, nine_zeros);
						}
					}
				}

				if (has_ctf_premultiplied && has_not_ctf_premultiplied)
				{
					if (!MDoptics[i].containsLabel(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED))
					{
						FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
						{
							MDoptics[i].setValue(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, false);
						}
					}
				}
			}

			// Now combine all optics tables into one
			obsModel.opticsMdt = MetaDataTable::combineMetaDataTables(MDoptics);
		}

		// Combine the particles tables
		MDout = MetaDataTable::combineMetaDataTables(MDsin);

		//Deactivate the group_name column
		MDout.deactivateLabel(EMDL_MLMODEL_GROUP_NO);

		if (fn_check != "")
		{
			EMDLabel label = EMDL::str2Label(fn_check);
			if (!MDout.containsLabel(label))
				REPORT_ERROR("ERROR: the output file does not contain the label to check for duplicates. Is it present in all input files?");

			/// Don't want to mess up original order, so make a MDsort with only that label...
			FileName fn_this, fn_prev = "";
			MetaDataTable MDsort;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDout)
			{
				MDout.getValue(label, fn_this);
				MDsort.addObject();
				MDsort.setValue(label, fn_this);
			}
			// sort on the label
			MDsort.newSort(label);
			long int nr_duplicates = 0;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsort)
			{
				MDsort.getValue(label, fn_this);
				if (fn_this == fn_prev)
				{
					nr_duplicates++;
					std::cerr << " WARNING: duplicate entry: " << fn_this << std::endl;
				}
				fn_prev = fn_this;
			}

			if (nr_duplicates > 0)
				std::cerr << " WARNING: Total number of duplicate "<< fn_check << " entries: " << nr_duplicates << std::endl;
		}

		write_check_ignore_optics(MDout, fn_out, MDin.getName());
		std::cout << " Written: " << fn_out << std::endl;
	}

	void split()
	{
		MetaDataTable MD;
		read_check_ignore_optics(MD, fn_in);

		// Randomise if neccesary
		if (do_random_order)
		{
			if (random_seed < 0)
				randomize_random_generator();
			else
				init_random_generator(random_seed);

			MD.randomiseOrder();
		}

		long int n_obj = MD.numberOfObjects();
		if (n_obj == 0)
		{
			REPORT_ERROR("ERROR: empty STAR file...");
		}

		if (nr_split < 0 && size_split < 0)
		{
			REPORT_ERROR("ERROR: nr_split and size_split are both zero. Set at least one of them to be positive.");
		}
		else if (nr_split < 0 && size_split > 0)
		{
			nr_split = CEIL(1. * n_obj / size_split);
		}
		else if (nr_split > 0 && size_split < 0)
		{
			size_split = CEIL(1. * n_obj / nr_split);
		}

		std::vector<MetaDataTable > MDouts;
		MDouts.resize(nr_split);

		long int n = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			int my_split = n / size_split;
			if (my_split < nr_split)
			{
				MDouts[my_split].addObject(MD.getObject(current_object));
			}
			else
			{
				break;
			}
			n++;
		}

		// Sjors 19jun2019: write out a star file with the output nodes
		MetaDataTable MDnodes;
		MDnodes.setName("output_nodes");
		FileName fnt0;
		fnt0 = integerToString(nr_split);
		for (int isplit = 0; isplit < nr_split; isplit ++)
		{
			FileName fnt = fn_out.insertBeforeExtension("_split"+integerToString(isplit+1));
			write_check_ignore_optics(MDouts[isplit], fnt, MD.getName());
			std::cout << " Written: " <<fnt << " with " << MDouts[isplit].numberOfObjects() << " objects." << std::endl;

			MDnodes.addObject();
			MDnodes.setValue(EMDL_PIPELINE_NODE_NAME, fnt);
			int type;
			if (MD.getName() == "micrographs")
			{
				type = NODE_MICS;
			}
			else if (MD.getName() == "movies")
			{
				type = NODE_MOVIES;
			}
			else
			{
				// just assume these are particles
				type = NODE_PART_DATA;
			}

			MDnodes.setValue(EMDL_PIPELINE_NODE_TYPE, type);
		}

		// write out the star file with the output nodes
		FileName mydir = fn_out.beforeLastOf("/");
		if (mydir == "") mydir = ".";
		MDnodes.write(mydir + "/" + RELION_OUTPUT_NODES);

	}

	
////
	void do_bild()
	{
		EMDLabel label1=EMDL_ORIENT_ROT, label2=EMDL_ORIENT_TILT, label3=EMDL_ORIENT_PSI;
		MetaDataTable MD;
		read_check_ignore_optics(MD, fn_in);
		int MAX_=30000;
		int i;
		RFLOAT ROT[MAX_],TILT[MAX_],PDF[MAX_];
		for(i=0;i<MAX_;i++){
			ROT[i]=0.0;
			TILT[i]=0.0;
			PDF[i]=0.0;
		}
		int C_PDF=0;
		bool first_iteration=true;
		int count_iter=-1;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			count_iter+=1;
			if(count_iter%int(ROUND(MD.numberOfObjects()/10))==0 && count_iter>1){
				std::cout<<"Finish line "<<count_iter<<", has "<<C_PDF<<" sample points."<<std::endl;
			}
			RFLOAT val_rot,val_tilt,val_psi,rotp,tiltp,psip;
			MD.getValue(label1, val_rot);
			MD.getValue(label2, val_tilt);
			MD.getValue(label3, val_psi);
			if(first_iteration){
				ROT[C_PDF]=val_rot;
				TILT[C_PDF]=val_tilt;
				PDF[C_PDF]=1.;
				C_PDF+=1;
			}
			else{
				Matrix1D<RFLOAT> v1(3),v2(3);
				Euler_angles2direction(val_rot, val_tilt,v1);
				RFLOAT dot_value,distance;
				bool is_within=false;
				for (i=0;i<C_PDF;i++){
					Euler_angles2direction(ROT[i], TILT[i],v2);
					dot_value=dotProduct(v1,v2);
					if(dot_value>1.0)
						dot_value=1.0;
					if(dot_value<-1.0)
						dot_value=-1.0;
					distance=RAD2DEG(acos(dot_value));
					if(distance<bild_angle_step){
						is_within=true;
						PDF[i]+=1.;
						break;
					}
				}
				if(!is_within){
					
					ROT[C_PDF]=val_rot;
					TILT[C_PDF]=val_tilt;
					PDF[C_PDF]=1.;
					C_PDF+=1;
				}
				if(C_PDF>=MAX_-1){
					std::cout<<"Warning, in do_bild(), more points than allowed. Enlarge the angle step."<<std::endl;
				}
			}
			first_iteration=false;
		}
		RFLOAT SUM=0.0,SUM_SQR=0.0,PDF_MAX=-99999999.0;
		for (i=0;i<C_PDF;i++){
			if(PDF[i]>PDF_MAX)
				PDF_MAX=PDF[i];
			SUM+=PDF[i];
			SUM_SQR+=PDF[i]*PDF[i];
		}
		SUM=SUM/C_PDF;
		SUM_SQR=SUM_SQR/C_PDF;
		RFLOAT PDF_MEAN=SUM,PDF_SIGMA=sqrt(SUM_SQR-SUM*SUM);
		std::cout<<"debug PDF_MEAN="<<PDF_MEAN<<", PDF_SIGMA="<<PDF_SIGMA<<", PDF_MAX="<<PDF_MAX<<std::endl;
		RFLOAT Rmax_frac = 0.3,width_frac = 0.5,offset=bild_radius,R=offset;
		std::ofstream tmp_string(fn_out);
		for (i=0;i<C_PDF;i++){
			RFLOAT colscale = (PDF[i] - PDF_MEAN) / PDF_SIGMA;
			if(colscale>5.0)
				colscale=5.0;
			if(colscale<-1.0)
				colscale=-1.0;
			colscale=colscale/6.0;
			colscale+=1.0/6.0;
			RFLOAT Rp = R + Rmax_frac * R * (PDF[i]) / PDF_MAX;
			Matrix1D<RFLOAT> AV(3);
			Euler_angles2direction(ROT[i], TILT[i],AV);
			RFLOAT width=bild_angle_step/2.0;
			if(width>0.5)
				width=0.5;
			if (fabs((R - Rp) * (AV(0))) > 0.01 ||  fabs((R - Rp) * (AV(1))) > 0.01 || fabs((R - Rp) * (AV(2))) > 0.01){
				tmp_string<<".color "<<(colscale)<<" 0 "<<(1. - colscale)<<std::endl;
				tmp_string<<".cylinder "<<(R*AV(0)+(offset))<<" "<<(R*AV(1)+(offset))<<" "<<(R*AV(2)+(offset))<<" "<<(Rp *AV(0)+ offset) <<" "<<(Rp *AV(1)+ offset)<<" "<<(Rp *AV(2)+ offset)<<" "<<(width)<<std::endl;
			}
		}
		tmp_string.close();
		std::cout<<"Finish do_bild. EXIT."<<std::endl;
	}
////

void do_rot_star()
	{
		EMDLabel label1=EMDL_ORIENT_ROT, label2=EMDL_ORIENT_TILT, label3=EMDL_ORIENT_PSI;

		MetaDataTable MD;
		read_check_ignore_optics(MD, fn_in);
		Matrix2D<RFLOAT> rotate_3D_L,rotate_3D_R;
		Euler_rotation3DMatrix(rOt, tIlt, pSi, rotate_3D_R);
		rotate_3D_R.resize(3, 3);
		rotate_3D_R.inv(rotate_3D_L);
		rotate_3D_L.resize(3, 3);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{

			RFLOAT val_rot,val_tilt,val_psi,rotp,tiltp,psip;
			MD.getValue(label1, val_rot);
			MD.getValue(label2, val_tilt);
			MD.getValue(label3, val_psi);
			Matrix2D<RFLOAT> inmap_L,inmap_R,inmap_mult;
			Euler_rotation3DMatrix(val_rot, val_tilt, val_psi, inmap_R);
			inmap_R.resize(3, 3);
			inmap_R.inv(inmap_L);
			inmap_L.resize(3, 3);
	
			inmap_mult=inmap_R*rotate_3D_L;
		
			Euler_matrix2angles(inmap_mult,rotp,tiltp,psip);
		
			MD.setValue(label1,rotp);
			MD.setValue(label2,tiltp);
			MD.setValue(label3,psip);
		}

		write_check_ignore_optics(MD, fn_out, MD.getName());
		std::cout << "do_RoT finished. Written: " << fn_out << std::endl;
	}
////
////
	void do_block_based()
	{
		EMDLabel label1=EMDL_ORIENT_ROT, label2=EMDL_ORIENT_TILT, label3=EMDL_ORIENT_PSI;
		EMDLabel label4=EMDL_ORIENT_ORIGIN_X_ANGSTROM, label5=EMDL_ORIENT_ORIGIN_Y_ANGSTROM, label6=EMDL_IMAGE_COORD_X, label7=EMDL_IMAGE_COORD_Y;
		EMDLabel label8=EMDL_CTF_DEFOCUSU, label9=EMDL_CTF_DEFOCUSV;
		EMDLabel lable_BBR1=EMDL_BBR_DELTA_Z,lable_BBR2=EMDL_BBR_PARTICLE_SERIAL_NUM;
		MetaDataTable MD,Mout;
	//	read_check_ignore_optics(MD, fn_in);
		do_ignore_optics = false;
		ObservationModel::loadSafely(fn_in, obsModel, MD, "discover", 1, false); // false means don't die upon failure
		if (obsModel.opticsMdt.numberOfObjects() == 0)
		{
			do_ignore_optics = true;
			std::cout << " + WARNING: reading input STAR file without optics groups ..." << std::endl;
			MD.read(fn_in);
		}
		SymList SL;
		SL.read_sym_file(fn_sym);
		if (SL.SymsNo() < 1)
		//	REPORT_ERROR("ERROR Nothing to do. Provide a point group with symmetry!");
			std::cout<<"SL.SymsNo() < 1. Check whether symmetry == C1."<<std::endl;
		RFLOAT pixel_size;
		int image_size;
	
		obsModel.opticsMdt.getValue(EMDL_IMAGE_PIXEL_SIZE, pixel_size);
		
		obsModel.opticsMdt.getValue(EMDL_IMAGE_SIZE, image_size);
		std::cout<<"apix="<<pixel_size<<" , ysize="<<image_size<<std::endl;
		Mout.clear();
		Matrix2D<RFLOAT> L(3,3), R(3,3); // A matrix from the list
		
		int BBR_particle_serial=0,dropped_particles=0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			BBR_particle_serial+=1;
			RFLOAT val_rot,val_tilt,val_psi,val_oriXA,val_oriYA,val_coordX,val_coordY,val_dfu,val_dfv;
			RFLOAT TMP_DFU,TMP_DFV;
			RFLOAT rotp=0.0,tiltp=0.0,psip=0.0;
			MD.getValue(label1, val_rot);
			MD.getValue(label2, val_tilt);
			MD.getValue(label3, val_psi);
			MD.getValue(label4, val_oriXA);
			MD.getValue(label5, val_oriYA);
			MD.getValue(label6, val_coordX);
			MD.getValue(label7, val_coordY);
			MD.getValue(label8, val_dfu);
			MD.getValue(label9, val_dfv);
			TMP_DFU=val_dfu;
			TMP_DFV=val_dfv;
			if(SL.SymsNo()>-1)////always run this. because the SymList does not contain unity itself.
			{
				Matrix2D<RFLOAT> MuLt;
				Euler_angles2matrix(val_rot, val_tilt, val_psi,MuLt);
				RFLOAT Xposi,Yposi,Zposi;
				Matrix1D<RFLOAT> IS(3);
				Matrix1D<RFLOAT> RoW;
				IS.resize(3);
				IS(0)=Posi(0)-(float)image_size/2.0;
				IS(1)=Posi(1)-(float)image_size/2.0;
				IS(2)=Posi(2)-(float)image_size/2.0;
				MuLt.getRow(0,RoW);
				Xposi=dotProduct(RoW,IS);
				MuLt.getRow(1,RoW);
				Yposi=dotProduct(RoW,IS);
				MuLt.getRow(2,RoW);
				Zposi=dotProduct(RoW,IS);
				
				if(ROUND_SHIFT){
					Xposi=ROUND(Xposi);
					Yposi=ROUND(Yposi);
				}
				Xposi*=-pixel_size;
				Yposi*=-pixel_size;
				Zposi*=pixel_size;
				Xposi+=val_oriXA;
				Yposi+=val_oriYA;
				RFLOAT inmicrogrpah_X=val_coordX-Xposi/pixel_size;
				RFLOAT inmicrogrpah_Y=val_coordY-Yposi/pixel_size;
				RFLOAT DZ=0.0;
				if(handerness){
					val_dfu=TMP_DFU+Zposi;
					val_dfv=TMP_DFV+Zposi;
					DZ=Zposi;
				}
				if(handerness_negative){
					val_dfu=TMP_DFU-Zposi;
					val_dfv=TMP_DFV-Zposi;
					DZ=-1.0*Zposi;
				}
				if(do_microgrpah_drop){
					if(inmicrogrpah_X<0 || inmicrogrpah_X >micrograph_Xsize || inmicrogrpah_Y<0 || inmicrogrpah_Y >micrograph_Ysize){
						dropped_particles+=1;
						continue;
					}
				}
				
				Mout.addObject();
				Mout.setObject(MD.getObject());
				Mout.setValue(EMDL_ORIENT_ROT, val_rot);
				Mout.setValue(EMDL_ORIENT_TILT, val_tilt);
				Mout.setValue(EMDL_ORIENT_PSI, val_psi);
				Mout.setValue(label4,Xposi);
				Mout.setValue(label5,Yposi);
				Mout.setValue(label8,val_dfu);
				Mout.setValue(label9,val_dfv);
				Mout.setValue(lable_BBR1,DZ);
				Mout.setValue(lable_BBR2,BBR_particle_serial);
			}
			if(SL.SymsNo()>0){
				
				for (int isym = 0; isym < SL.SymsNo(); isym++)
				{
					
					SL.get_matrices(isym, L, R);
					L.resize(3, 3); // Erase last row and column
					R.resize(3, 3); // as only the relative orientation is useful and not the translation

					Euler_apply_transf(L, R, val_rot, val_tilt, val_psi, rotp, tiltp, psip);
					
					Matrix2D<RFLOAT> NEW_INMAP,INMAP,INV_INMAP,MuLt;
					Matrix1D<RFLOAT> RoW;
					Euler_angles2matrix(rotp, tiltp, psip,NEW_INMAP);
					RFLOAT Xposi,Yposi,Zposi;
					MuLt=NEW_INMAP;
					Matrix1D<RFLOAT> IS(3);
					IS.resize(3);
					IS(0)=Posi(0)-(float)image_size/2.0;
					IS(1)=Posi(1)-(float)image_size/2.0;
					IS(2)=Posi(2)-(float)image_size/2.0;
					MuLt.getRow(0,RoW);
					Xposi=dotProduct(RoW,IS);
					MuLt.getRow(1,RoW);
					Yposi=dotProduct(RoW,IS);
					MuLt.getRow(2,RoW);
					Zposi=dotProduct(RoW,IS);
					if(ROUND_SHIFT){
						Xposi=ROUND(Xposi);
						Yposi=ROUND(Yposi);
					}
					Xposi*=-pixel_size;
					Yposi*=-pixel_size;
					Zposi*=pixel_size;
					Xposi+=val_oriXA;
					Yposi+=val_oriYA;
					RFLOAT inmicrogrpah_X=val_coordX-Xposi/pixel_size;
					RFLOAT inmicrogrpah_Y=val_coordY-Yposi/pixel_size;
					RFLOAT DZ=0.0;
					if(handerness){
						val_dfu=TMP_DFU+Zposi;
						val_dfv=TMP_DFV+Zposi;
						DZ=Zposi;
					}
					if(handerness_negative){
						val_dfu=TMP_DFU-Zposi;
						val_dfv=TMP_DFV-Zposi;
						DZ=-1.0*Zposi;
					}
					if(do_microgrpah_drop){
						if(inmicrogrpah_X<0 || inmicrogrpah_X >micrograph_Xsize || inmicrogrpah_Y<0 || inmicrogrpah_Y >micrograph_Ysize){
							dropped_particles+=1;
							continue;
						}
					}
					
					Mout.addObject();
					Mout.setObject(MD.getObject());
					Mout.setValue(EMDL_ORIENT_ROT, rotp);
					Mout.setValue(EMDL_ORIENT_TILT, tiltp);
					Mout.setValue(EMDL_ORIENT_PSI, psip);
					Mout.setValue(label4,Xposi);
					Mout.setValue(label5,Yposi);
					Mout.setValue(label8,val_dfu);
					Mout.setValue(label9,val_dfv);
					Mout.setValue(lable_BBR1,DZ);
					Mout.setValue(lable_BBR2,BBR_particle_serial);
				
				}
			}
		}
		write_check_ignore_optics(Mout, fn_out, MD.getName());
		std::cout << "Block calculation finished. Written: " << fn_out << std::endl;
		if(do_microgrpah_drop){
			std::cout << "Total number of " << dropped_particles <<" particles were dropped dueing to out of the micrograph."<< std::endl;
		}

	}
////
	void invert_BBR_handerness()
	{
		EMDLabel lable_BBR1=EMDL_BBR_DELTA_Z;
		EMDLabel label1=EMDL_CTF_DEFOCUSU, label2=EMDL_CTF_DEFOCUSV;
		MetaDataTable MD;
		read_check_ignore_optics(MD, fn_in);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			RFLOAT val_dz,val_dfu,val_dfv,new_dfu,new_dfv;
			MD.getValue(label1, val_dfu);
			MD.getValue(label2, val_dfv);
			MD.getValue(lable_BBR1, val_dz);
			new_dfu=val_dfu-2.0*val_dz;
			new_dfv=val_dfv-2.0*val_dz;
			MD.setValue(label1,new_dfu);
			MD.setValue(label2,new_dfv);
		}

		write_check_ignore_optics(MD, fn_out, MD.getName());
		std::cout << "do_invert_BBR_handerness finished. Written: " << fn_out << std::endl;
	}
////
	void operate()
	{
		EMDLabel label1, label2, label3;
		label1 = EMDL::str2Label(fn_operate);
		if (fn_operate2 != "")
		{
			label2 = EMDL::str2Label(fn_operate2);
		}
		if (fn_operate3 != "")
		{
			label3 = EMDL::str2Label(fn_operate3);
		}

		MetaDataTable MD;
		read_check_ignore_optics(MD, fn_in);

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			if (EMDL::isDouble(label1))
			{
				RFLOAT val;
				if (fn_set != "")
				{
					val = textToFloat(fn_set);
					MD.setValue(label1, val);
					if (fn_operate2 != "") MD.setValue(label2, val);
					if (fn_operate3 != "") MD.setValue(label3, val);
				}
				else if (multiply_by != 1. || add_to != 0.)
				{
					MD.getValue(label1, val);
					val = multiply_by * val + add_to;
					MD.setValue(label1, val);
					if (fn_operate2 != "")
					{
						MD.getValue(label2, val);
						val = multiply_by * val + add_to;
						MD.setValue(label2, val);
					}
					if (fn_operate3 != "")
					{
						MD.getValue(label3, val);
						val = multiply_by * val + add_to;
						MD.setValue(label3, val);
					}
				}
			}
			else if (EMDL::isInt(label1))
			{
				int val;
				if (fn_set != "")
				{
					val = textToInteger(fn_set);
					MD.setValue(label1, val);
					if (fn_operate2 != "") MD.setValue(label2, val);
					if (fn_operate3 != "") MD.setValue(label3, val);
				}
				else if (multiply_by != 1. || add_to != 0.)
				{
					MD.getValue(label1, val);
					val = multiply_by * val + add_to;
					MD.setValue(label1, val);
					if (fn_operate2 != "")
					{
						MD.getValue(label2, val);
						val = multiply_by * val + add_to;
						MD.setValue(label2, val);
					}
					if (fn_operate3 != "")
					{
						MD.getValue(label3, val);
						val = multiply_by * val + add_to;
						MD.setValue(label3, val);
					}
				}

			}
			else if (EMDL::isString(label1))
			{
				if (fn_set != "")
				{
					MD.setValue(label1, fn_set);
					if (fn_operate2 != "") MD.setValue(label2, fn_set);
					if (fn_operate3 != "") MD.setValue(label3, fn_set);
				}
				else if (multiply_by != 1. || add_to != 0.)
				{
					REPORT_ERROR("ERROR: cannot multiply_by or add_to a string!");
				}
			}
			else if (EMDL::isBool(label1))
			{
				REPORT_ERROR("ERROR: cannot operate on a boolean!");
			}
			else if (EMDL::isBool(label1))
			{
				// @TODO:
				REPORT_ERROR("ERROR: cannot operate on vectors (yet)!");
			}

		}

		write_check_ignore_optics(MD, fn_out, MD.getName());
		std::cout << " Written: " << fn_out << std::endl;
	}

	void center()
	{
		MetaDataTable MD;
		read_check_ignore_optics(MD, fn_in, "particles");
		bool do_contains_xy = (MD.containsLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM) && MD.containsLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM));
		bool do_contains_z = (MD.containsLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM));

		if (!do_contains_xy)
		{
			REPORT_ERROR("ERROR: input STAR file does not contain rlnOriginX/Y for re-centering.");
		}

		Matrix1D<RFLOAT> my_center(3);
		XX(my_center) = center_X;
		YY(my_center) = center_Y;
		ZZ(my_center) = center_Z;

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			int optics_group;
			Matrix1D<RFLOAT> my_projected_center(3);
			Matrix2D<RFLOAT> A3D;
			RFLOAT xoff, yoff, zoff, rot, tilt, psi, angpix;

			if (do_ignore_optics)
			{
				angpix = cl_angpix;
			}
			else
			{
				MD.getValue(EMDL_IMAGE_OPTICS_GROUP, optics_group);
				optics_group--;
				angpix = obsModel.getPixelSize(optics_group);
			}

			MD.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff);
			MD.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff);
			MD.getValue(EMDL_ORIENT_ROT, rot);
			MD.getValue(EMDL_ORIENT_TILT, tilt);
			MD.getValue(EMDL_ORIENT_PSI, psi);

			xoff /= angpix;
			yoff /= angpix;

			// Project the center-coordinates
			Euler_angles2matrix(rot, tilt, psi, A3D, false);
			my_projected_center = A3D * my_center;

			xoff -= XX(my_projected_center);
			yoff -= YY(my_projected_center);

			// Set back the new centers
			MD.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff*angpix);
			MD.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff*angpix);

			// also allow 3D data (subtomograms)
			if (do_contains_z)
			{
				MD.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zoff);
				zoff /= angpix;
				zoff -= ZZ(my_projected_center);
				MD.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zoff*angpix);
			}
		}

		write_check_ignore_optics(MD, fn_out, MD.getName());
		std::cout << " Written: " << fn_out << std::endl;
	}

	void remove_column()
	{
		MetaDataTable MD;
		read_check_ignore_optics(MD, fn_in);
		MD.deactivateLabel(EMDL::str2Label(remove_col_label));
		write_check_ignore_optics(MD, fn_out, MD.getName());
		std::cout << " Written: " << fn_out << std::endl;
	}

	void add_column()
	{
		if ((add_col_value == "" && add_col_from == "") ||
		    (add_col_value != "" && add_col_from != ""))
			REPORT_ERROR("ERROR: you need to specify either --add_column_value or --copy_column_from when adding a column.");

		bool set_value = (add_col_value != "");

		MetaDataTable MD;
		EMDLabel label = EMDL::str2Label(add_col_label);
		EMDLabel source_label;

		read_check_ignore_optics(MD, fn_in);
		MD.addLabel(label);

		if (add_col_from != "")
		{
			source_label = EMDL::str2Label(add_col_from);
			if (!MD.containsLabel(source_label))
				REPORT_ERROR("ERROR: The column specified in --add_column_from is not present in the input STAR file.");
		}

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			if (EMDL::isDouble(label))
			{
				RFLOAT aux;
				if (set_value) aux = textToFloat(add_col_value);
				else MD.getValue(source_label, aux);

				MD.setValue(label, aux);
			}
			else if (EMDL::isInt(label))
			{
				long aux;
				if (set_value) aux = textToInteger(add_col_value);
				else MD.getValue(source_label, aux);

				MD.setValue(label, aux);
			}
			else if (EMDL::isBool(label))
			{
				bool aux;
				if (set_value) aux = (bool)textToInteger(add_col_value);
				else MD.getValue(source_label, aux);

				MD.setValue(label, aux);
			}
			else if (EMDL::isString(label))
			{
				std::string aux;
				if (set_value) aux = add_col_value;
				else MD.getValue(source_label, aux);

				MD.setValue(label, add_col_value);
			}
			else if (EMDL::isString(label))
			{
				std::string auxStr;
				if (set_value) auxStr = add_col_value;
				else MD.getValueToString(source_label, auxStr);

				MD.setValueFromString(label, add_col_value);
			}
		}

		write_check_ignore_optics(MD, fn_out, MD.getName());
		std::cout << " Written: " << fn_out << std::endl;
	}

	void hist_column()
	{
		MetaDataTable MD;
		EMDLabel label = EMDL::str2Label(hist_col_label);

		std::vector<RFLOAT> values;

		read_check_ignore_optics(MD, fn_in);
		if (!MD.containsLabel(label))
			REPORT_ERROR("ERROR: The column specified in --hist_column is not present in the input STAR file.");

		std::vector<RFLOAT> histX,histY;
		CPlot2D *plot2D=new CPlot2D("");
		MD.columnHistogram(label, histY, histX, 1, plot2D, nr_bin, hist_min, hist_max, show_frac, show_cumulative);
		FileName fn_eps = fn_out.withoutExtension()+".eps";
		plot2D->OutputPostScriptPlot(fn_eps);
		std::cout << " Done! written out " << fn_eps << std::endl;
		delete plot2D;

	}

	void remove_duplicate()
	{
		if (do_ignore_optics)
			REPORT_ERROR("Duplicate removal is not compatible with --ignore_optics");

		MetaDataTable MD;
		read_check_ignore_optics(MD, fn_in, "particles");

		EMDLabel mic_label;
		if (MD.containsLabel(EMDL_MICROGRAPH_NAME)) mic_label = EMDL_MICROGRAPH_NAME;
		else REPORT_ERROR("The input STAR file does not contain rlnMicrographName column.");

		RFLOAT particle_angpix = 1.0; // rlnOriginX/YAngst is always 1 A/px.

		if (obsModel.numberOfOpticsGroups() > 1)
			std::cerr << "WARNING: The input contains multiple optics groups. We assume that the pixel sizes of original micrographs before extraction are all the same. If this is not the case, you have to split the input and remove duplicates separately." << std::endl;

		if (extract_angpix > 0)
		{
			std::cout << " + Using the provided pixel size for original micrographs before extraction: " << extract_angpix << std::endl;
		}
		else
		{
			extract_angpix = obsModel.getPixelSize(0);
			std::cout << " + Assuming the pixel size of original micrographs before extraction is " << extract_angpix << std::endl;
		}

		RFLOAT scale = particle_angpix / extract_angpix;
		RFLOAT duplicate_threshold_in_px = duplicate_threshold / extract_angpix;

		std::cout << " + The minimum inter-particle distance " << duplicate_threshold << " A corresponds to " << duplicate_threshold_in_px << " px in the micrograph coordinate (rlnCoordinateX/Y)." << std::endl;
		std::cout << " + The particle shifts (rlnOriginXAngst, rlnOriginYAngst) are multiplied by " << scale << " to bring it to the same scale as rlnCoordinateX/Y." << std::endl;
		FileName fn_removed = fn_out.withoutExtension() + "_removed.star";


		MetaDataTable MDout = removeDuplicatedParticles(MD, mic_label, duplicate_threshold_in_px, scale, fn_removed, true);

		write_check_ignore_optics(MDout, fn_out, "particles");
		std::cout << " Written: " << fn_out << std::endl;
	}
};

int main(int argc, char *argv[])
{
	star_handler_parameters prm;

	try
	{
		prm.read(argc, argv);
		prm.run();

	}
	catch (RelionError XE)
	{
        	std::cerr << XE;
		//prm.usage();
		return RELION_EXIT_FAILURE;
	}
	return RELION_EXIT_SUCCESS;
}

