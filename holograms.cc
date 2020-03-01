//-----Richard Luis Martin 2012/12/15
//-----This code generates so-called Voronoi holograms, which are 3D vector
//-----representations of the accessible void space within a material.
//-----With respect to a certain probe size, the accessible voronoi network
//-----is calculated and each edge of the network is characterized in terms of
//-----length, start and end radii. The hologram stores the count of each edge
//-----type with respect to these three properties.
//-----The code outputs: 1) "*.stats", which describes the contents of the
// hologram; 2) "*_holo.txt" which contains only the non-zero entries; 3)
//"*_complete_holo.txt" which contains the complete holograms written on one
// line.

//#include "network.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "holograms.h"

using namespace std;

// the main code for analyzing the accessible network and categorizing the
// edges, producing hologram representations
void analyze_accessible_voronoi_pre_segment(VORONOI_NETWORK *vornet,
                                            float probeRad,
                                            vector<bool> *accessInfo,
                                            char *name,
                                            const char *bin_directory) {

  // 1) inspect the accessibility of the network - specifically, we need to find
  // the accessible edges so they can be categorized

  cout << "NUM_BINS = " << NUM_BINS << "\n";

  // find info and count accessible nodes
  unsigned int i, j, k;
  unsigned int node_count = 0, edge_count = 0;
  float min_node_radius = -1, max_node_radius = -1, average_node_radius = 0,
        min_edge_radius = -1, max_edge_radius = -1, average_edge_radius = 0,
        min_edge_length = -1, max_edge_length = -1, average_edge_length = 0;
  vector<VOR_NODE> accessible_voronoi_nodes = vector<VOR_NODE>();
  for (i = 0; i < accessInfo->size(); i++) {
    if (accessInfo->at(i) || IS_FULL_VORONOI == 1) {
      average_node_radius += vornet->nodes.at(i).rad_stat_sphere;
      accessible_voronoi_nodes.push_back(vornet->nodes.at(i));
      if (vornet->nodes.at(i).rad_stat_sphere < min_node_radius ||
          min_node_radius < 0)
        min_node_radius = vornet->nodes.at(i).rad_stat_sphere;
      if (vornet->nodes.at(i).rad_stat_sphere > max_node_radius ||
          max_node_radius < 0)
        max_node_radius = vornet->nodes.at(i).rad_stat_sphere;
      node_count++;
    }
  }
  cout << node_count << " accessible nodes in total."
       << "\n";
  if (node_count > 0) {
    average_node_radius /= node_count;
  }
  // based on accessible nodes, determine accessible edges
  vector<VOR_EDGE> accessible_voronoi_edges = vector<VOR_EDGE>();
  vector<unsigned int> accessible_voronoi_edge_indices = vector<unsigned int>();
  for (i = 0; i < vornet->edges.size(); i++) {
    if ((accessInfo->at(vornet->edges.at(i).from) &&
         accessInfo->at(vornet->edges.at(i).to) &&
         vornet->edges.at(i).rad_moving_sphere > probeRad) ||
        IS_FULL_VORONOI ==
            1) { // i.e. if both nodes for this edge are accessible, and if the
                 // radius is larger than the probe radius
      average_edge_radius += vornet->edges.at(i).rad_moving_sphere;
      average_edge_length += vornet->edges.at(i).length;
      if (vornet->edges.at(i).rad_moving_sphere < min_edge_radius ||
          min_edge_radius < 0)
        min_edge_radius = vornet->edges.at(i).rad_moving_sphere;
      if (vornet->edges.at(i).rad_moving_sphere > max_edge_radius ||
          max_edge_radius < 0)
        max_edge_radius = vornet->edges.at(i).rad_moving_sphere;
      if (vornet->edges.at(i).length < min_edge_length || min_edge_length < 0)
        min_edge_length = vornet->edges.at(i).length;
      if (vornet->edges.at(i).length > max_edge_length || max_edge_length < 0)
        max_edge_length = vornet->edges.at(i).length;
      accessible_voronoi_edges.push_back(vornet->edges.at(i));
      accessible_voronoi_edge_indices.push_back(i);
      edge_count++;
    }
  }
  cout << edge_count << " accessible edges in total."
       << "\n";
  if (edge_count > 0) {
    average_edge_radius /= edge_count;
    average_edge_length /= edge_count;
  }

  // write stats to terminal
  cout << "min_node_radius = " << min_node_radius
       << " max_node_radius = " << max_node_radius
       << " average_node_radius = " << average_node_radius << "\n";
  //	cout << "min_edge_radius = " << min_edge_radius << " max_edge_radius = "
  //<< max_edge_radius << " average_edge_radius = " << average_edge_radius <<
  //"\n";
  cout << "min_edge_length = " << min_edge_length
       << " max_edge_length = " << max_edge_length
       << " average_edge_length = " << average_edge_length << "\n";

  // 2) now we have the accessible edges, categorize them in terms of length,
  // start and end radii, using three separate histograms

  // now bin the nodes and edges following a pre-set binning system
  int int_node_radii_bins[NUM_BINS], int_edge_radii_bins[NUM_BINS],
      int_edge_length_bins[NUM_BINS];
  float float_node_radii_bins[NUM_BINS], float_edge_radii_bins[NUM_BINS],
      float_edge_length_bins[NUM_BINS];
  // geometric progression custom binning for edge lengths
  float *float_node_radii_bin_upper_bounds, *float_edge_radii_bin_upper_bounds,
      *float_edge_length_bin_upper_bounds;
  float_node_radii_bin_upper_bounds = new float[NUM_BINS],
  float_edge_radii_bin_upper_bounds = new float[NUM_BINS],
  float_edge_length_bin_upper_bounds = new float[NUM_BINS];
  // a 3D binning system considering for each edge its 1) length, 2/3) radii of
  // largest/smallest attached node
  int int_edge_stats_bins[NUM_BINS][NUM_BINS][NUM_BINS];
  float float_edge_stats_bins[NUM_BINS][NUM_BINS][NUM_BINS];

  // BIN READING
  // this process does one of two things:
  // 1) no bin_directory is specified - the default bins are used
  // 2) a bin_directory is given - grab the files from there
  if (!bin_directory) { // case 1)
    printf("Bin directory not specified - using defaults.\n");
    for (i = 0; i < NUM_BINS - 1; i++) {
      float_edge_length_bin_upper_bounds[i] = default_edge_length_bins[i];
      float_node_radii_bin_upper_bounds[i] = default_node_radii_bins[i];
    }
  } else { // case 2)
    FILE *edge_length_bin_file, *node_radii_bin_file;
    char *edge_length_bin_string = new char[100];
    strcpy(edge_length_bin_string, bin_directory);
    strcat(edge_length_bin_string, "edge_length_bins");
    char *node_radii_bin_string = new char[100];
    strcpy(node_radii_bin_string, bin_directory);
    strcat(node_radii_bin_string, "node_radii_bins");
    printf("Bin directory specified - trying to open files %s and %s.\n",
           edge_length_bin_string, node_radii_bin_string);
    edge_length_bin_file = fopen(edge_length_bin_string, "r");
    node_radii_bin_file = fopen(node_radii_bin_string, "r");
    char error = 0;
    if (edge_length_bin_file == NULL) {
      printf("ERROR: could not open bins file %s\n", edge_length_bin_string);
      error = 1;
    }
    if (node_radii_bin_file == NULL) {
      printf("ERROR: could not open bins file %s\n", node_radii_bin_string);
      error = 1;
    }
    if (error == 1)
      exit(EXIT_FAILURE);
    printf("Files for bin bounds opened successfully.\n");
    // read in the values for each bin from the files
    for (i = 0; i < NUM_BINS - 1; i++) {
      float temp1, temp2;
      int status;
      status = fscanf(edge_length_bin_file, "%f", &temp1);
      float_edge_length_bin_upper_bounds[i] = temp1;
      if (status == -1) {
        printf("ERROR: could not read edge length bin bound number %d from "
               "file %s\n",
               i + 1, edge_length_bin_string);
        exit(EXIT_FAILURE);
      }
      status = fscanf(node_radii_bin_file, "%f", &temp2);
      float_node_radii_bin_upper_bounds[i] = temp2;
      if (status == -1) {
        printf("ERROR: could not read node radii bin bound number %d from file "
               "%s\n",
               i + 1, node_radii_bin_string);
        exit(EXIT_FAILURE);
      }
    }
    fclose(edge_length_bin_file);
    fclose(node_radii_bin_file);
    // we've read the bins in so can make the holograms
  }
  // now check that the bins make sense
  int error = 0;
  for (i = 1; i < NUM_BINS - 1 && error == 0; i++) {
    if (float_edge_length_bin_upper_bounds[i] <
        float_edge_length_bin_upper_bounds[i - 1])
      error = 1;
    if (error == 1) {
      printf("ERROR: edge length bins do not constitute a non-decreasing "
             "series (entry index %d = %f; entry index %d = %f)\n",
             i - 1, float_edge_length_bin_upper_bounds[i - 1], i,
             float_edge_length_bin_upper_bounds[i]);
      if (!bin_directory)
        printf("\t\tPlease contact the developers since this problem occurred "
               "for the default bins\n");
      exit(EXIT_FAILURE);
    }
    if (float_node_radii_bin_upper_bounds[i] <
        float_node_radii_bin_upper_bounds[i - 1])
      error = 1;
    if (error == 1) {
      printf("ERROR: node radii bins do not constitute a non-decreasing series "
             "(entry index %d = %f; entry index %d = %f)\n",
             i - 1, float_node_radii_bin_upper_bounds[i - 1], i,
             float_node_radii_bin_upper_bounds[i]);
      if (!bin_directory)
        printf("\t\tPlease contact the developers since this problem occurred "
               "for the default bins\n");
      exit(EXIT_FAILURE);
    }
  }

  // 3) prepare the 3D vectors that will store the hologram data, and fill them
  // now safe to move on
  float_node_radii_bin_upper_bounds[NUM_BINS - 1] = 1000;  // i.e. 'infinity'
  float_edge_length_bin_upper_bounds[NUM_BINS - 1] = 1000; // i.e. 'infinity'
  for (i = 0; i < NUM_BINS; i++) {
    int_node_radii_bins[i] = 0;
    int_edge_length_bins[i] = 0;
    float_node_radii_bins[i] = 0;
    float_edge_length_bins[i] = 0;
    for (j = 0; j < NUM_BINS; j++) {
      for (k = 0; k < NUM_BINS; k++) {
        int_edge_stats_bins[i][j][k] = 0;
        float_edge_stats_bins[i][j][k] = 0;
      }
    }
  }
  // only do anything if there are any nodes!
  if (node_count > 0) {
    for (i = 0; i < accessible_voronoi_nodes.size(); i++) {
      float radius = accessible_voronoi_nodes.at(i).rad_stat_sphere;
      int bin_index = -1;
      bin_index = get_bin(radius, float_node_radii_bin_upper_bounds);
      if (IS_COUNT == 1)
        int_node_radii_bins[bin_index]++; // 1 for count
      else
        int_node_radii_bins[bin_index] += 100; // 100 so we can get percentages
    }
    for (i = 0; i < accessible_voronoi_edges.size(); i++) {
      float radius = accessible_voronoi_edges.at(i).rad_moving_sphere,
            length = accessible_voronoi_edges.at(i).length;
      int bin_index = -1;
      bin_index = get_bin(radius, float_edge_radii_bin_upper_bounds);
      if (IS_COUNT == 1)
        int_edge_radii_bins[bin_index]++; // 1 for count
      else
        int_edge_radii_bins[bin_index] += 100; // 100 so we can get percentages
      bin_index = -1;
      bin_index = get_bin(length, float_edge_length_bin_upper_bounds);
      if (IS_COUNT == 1)
        int_edge_length_bins[bin_index]++; // 1 for count
      else
        int_edge_length_bins[bin_index] += 100; // 100 so we can get percentages
      // now handle the 3D edge-based descriptors
      float from_radii =
          vornet->nodes
              .at(vornet->edges.at(accessible_voronoi_edge_indices.at(i)).from)
              .rad_stat_sphere;
      float to_radii =
          vornet->nodes
              .at(vornet->edges.at(accessible_voronoi_edge_indices.at(i)).to)
              .rad_stat_sphere;
      float min_radii = min(from_radii, to_radii),
            max_radii = max(from_radii, to_radii);

      if (IS_COUNT == 1)
        int_edge_stats_bins
            [get_bin(length, float_edge_length_bin_upper_bounds)]
            [get_bin(max_radii, float_node_radii_bin_upper_bounds)][get_bin(
                min_radii, float_node_radii_bin_upper_bounds)]++; // 1 for count
      else
        int_edge_stats_bins[get_bin(length, float_edge_length_bin_upper_bounds)]
                           [get_bin(max_radii,
                                    float_node_radii_bin_upper_bounds)]
                           [get_bin(min_radii,
                                    float_node_radii_bin_upper_bounds)] +=
            100; // 100 so we can get percentages
    }
    // now convert binning sums into binning proportions
    for (i = 0; i < NUM_BINS; i++) {
      float_node_radii_bins[i] =
          ((float)int_node_radii_bins[i]) / accessible_voronoi_nodes.size();
      float_edge_radii_bins[i] =
          ((float)int_edge_radii_bins[i]) / accessible_voronoi_edges.size();
      float_edge_length_bins[i] =
          ((float)int_edge_length_bins[i]) / accessible_voronoi_edges.size();
      for (j = 0; j < NUM_BINS; j++) {
        for (k = 0; k <= j; k++) { // because k is never larger than j
          float_edge_stats_bins[i][j][k] =
              ((float)int_edge_stats_bins[i][j][k]) /
              accessible_voronoi_edges.size();
        }
      }
    }
  } // end if any nodes

  // write stats to terminal - at this stage if no nodes, all percentages will
  // be zero
  if (IS_COUNT == 1) {
    cout << "edge_length bin counts:";
    for (i = 0; i < NUM_BINS; i++) {
      cout << " " << int_edge_length_bins[i];
    }
    cout << "\n";
    cout << "node_radius bin counts:";
    for (i = 0; i < NUM_BINS; i++) {
      cout << " " << int_node_radii_bins[i];
    }
    cout << "\n";
    //		cout << "edge_radius bin counts:";
    //		for(i=0; i<NUM_BINS; i++) {
    //		        cout << " " << int_edge_radii_bins[i];
    //		}
    //		cout << "\n";
  } else {
    cout << "node_radius bin percentages:";
    for (i = 0; i < NUM_BINS; i++) {
      cout << " " << float_node_radii_bins[i];
    }
    cout << "\n";
    //		cout << "edge_radius bin percentages:";
    //		for(i=0; i<NUM_BINS; i++) {
    //			cout << " " << float_edge_radii_bins[i];
    //		}
    //		cout << "\n";
    cout << "edge_length bin percentages:";
    for (i = 0; i < NUM_BINS; i++) {
      cout << " " << float_edge_length_bins[i];
    }
    cout << "\n";
  }

  // 4) file output section; writes three files: "*.stats", which describes the
  // contents of the hologram; "*.holo.txt", which is a concise representation
  // storing only the non-zero entries, and also in a format which can be
  // visualized in the VisIt software package; "*complete_holo.txt", which is a
  // one-line print out of the complete 3D hologram, including zero-valued
  // entries.

  // now write to output file
  char outputFile[256];
  strcpy(outputFile, name);
  strcat(outputFile, ".stats");
  FILE *output;
  output = fopen(outputFile, "w");
  if (output == NULL) {
    printf("ERROR: cannot open file %s to write Voronoi hologram statistics\n",
           outputFile);
    exit(EXIT_FAILURE);
  }
  fprintf(output,
          "%s statistics:\n %d Nodes:\n  min_radius %.6f\n  max_radius %.6f\n  "
          "average_radius %.6f\n %d Edges:\n  min_radius %.6f\n  max_radius "
          "%.6f\n  average_radius %.6f\n  min_length %.6f\n  max_length %.6f\n "
          " average_length %.6f\n",
          name, node_count, min_node_radius, max_node_radius,
          average_node_radius, edge_count, min_edge_radius, max_edge_radius,
          average_edge_radius, min_edge_length, max_edge_length,
          average_edge_length);
  if (IS_COUNT == 1) {
    fprintf(output, " Node_radii_bin_values:\n ");
    for (i = 0; i < NUM_BINS; i++) {
      fprintf(output, " %d", int_node_radii_bins[i]);
    }
    fprintf(output, "\n");
    fprintf(output, " Edge_length_bin_values:\n ");
    for (i = 0; i < NUM_BINS; i++) {
      fprintf(output, " %d", int_edge_length_bins[i]);
    }
    fprintf(output, "\n");
    fclose(output);
  } else {
    fprintf(output, " Node_radii_bin_values:\n ");
    for (i = 0; i < NUM_BINS; i++) {
      fprintf(output, " %.3f", float_node_radii_bins[i]);
    }
    fprintf(output, "\n");
    fprintf(output, " Edge_length_bin_values:\n ");
    for (i = 0; i < NUM_BINS; i++) {
      fprintf(output, " %.3f", float_edge_length_bins[i]);
    }
    fprintf(output, "\n");
    fclose(output);
  }
  /*
          fstream output; output.open(outputFile, fstream::out);
          output << name << ": " << node_count << " " << min_node_radius << " "
     << max_node_radius << " " << average_node_radius << " " << edge_count << "
     " << min_edge_radius << " " << max_edge_radius << " " <<
     average_edge_radius << " " << min_edge_length << " " << max_edge_length <<
     " " << average_edge_length << " "; if(IS_COUNT==1) { for(i=0; i<NUM_BINS;
     i++) { output << " " << int_node_radii_bins[i];
                  }
                  output << " ";
                  for(i=0; i<NUM_BINS; i++) {
                          output << " " << int_edge_radii_bins[i];
                  }
                  output << " ";
                  for(i=0; i<NUM_BINS; i++) {
                          output << " " << int_edge_length_bins[i];
                  }
                  output << "\n";
                  output.close();
          } else {
                  for(i=0; i<NUM_BINS; i++) {
                          output << " " << float_node_radii_bins[i];
                  }
                  output << " ";
                  for(i=0; i<NUM_BINS; i++) {
                          output << " " << float_edge_radii_bins[i];
                  }
                  output << " ";
                  for(i=0; i<NUM_BINS; i++) {
                          output << " " << float_edge_length_bins[i];
                  }
                  output << "\n";
                  output.close();
          }
  */

  // now write to second output file!
  char outputFile2[256];
  strcpy(outputFile2, name);
  strcat(outputFile2, "_holo.txt");
  /*  int num_bins_set = 0;
          for(i=0; i<NUM_BINS; i++) {
                  for(j=0; j<NUM_BINS; j++) {
                          for(k=0; k<=j; k++) { //because k is never larger than
    j if(int_edge_stats_bins[i][j][k] > 0) num_bins_set++;
        }
      }
    }*/
  FILE *output2;
  output2 = fopen(outputFile2, "w");
  if (output2 == NULL) {
    printf("ERROR: cannot open file %s to write Voronoi hologram xyz format "
           "file\n",
           outputFile2);
    exit(EXIT_FAILURE);
  }
  //	fstream output2; output2.open(outputFile2, fstream::out);
  //	output2 << num_bins_set << "\n" << name << " Voronoi hologram as xyz
  // format\n";
  fprintf(output2, "x y z frequency\n");
  if (IS_COUNT == 1) {
    for (i = 0; i < NUM_BINS; i++) {
      for (j = 0; j < NUM_BINS; j++) {
        for (k = 0; k <= j; k++) { // because k is never larger than j
          if (int_edge_stats_bins[i][j][k] > 0) {
            //            output2 << "# " << i << " " << j << " " << k << " " <<
            //            int_edge_stats_bins[i][j][k] << "\n";
            fprintf(output2, "%d %d %d %d\n", i, j, k,
                    int_edge_stats_bins[i][j][k]);
          }
        }
      }
    }
    fprintf(output2, "\n");
    fclose(output2);
  } else {
    for (i = 0; i < NUM_BINS; i++) {
      for (j = 0; j < NUM_BINS; j++) {
        for (k = 0; k <= j; k++) { // because k is never larger than j
          if (int_edge_stats_bins[i][j][k] > 0) {
            //						output2 << "# " << i << " " << j << " " << k << " "
            //<< float_edge_stats_bins[i][j][k] << "\n";
            fprintf(output2, "%d %d %d %.3f\n", i, j, k,
                    float_edge_stats_bins[i][j][k]);
          }
        }
      }
    }
    fprintf(output2, "\n");
    fclose(output2);
  }

  // now write to third output file!
  char outputFile3[256];
  strcpy(outputFile3, name);
  strcat(outputFile3, "_complete_holo.txt");
  FILE *output3;
  output3 = fopen(outputFile3, "w");
  if (output3 == NULL) {
    printf("ERROR: cannot open file %s to write Voronoi hologram complete data "
           "file\n",
           outputFile3);
    exit(EXIT_FAILURE);
  }
  //	fstream output3; output3.open(outputFile3, fstream::out);
  fprintf(output3, "%s:", name);
  if (IS_COUNT == 1) {
    for (i = 0; i < NUM_BINS; i++) {
      for (j = 0; j < NUM_BINS; j++) {
        for (k = 0; k <= j; k++) { // because k is never larger than j
          fprintf(output3, " %d", int_edge_stats_bins[i][j][k]);
        }
      }
    }
    fprintf(output3, "\n");
    fclose(output3);
  } else {
    for (i = 0; i < NUM_BINS; i++) {
      for (j = 0; j < NUM_BINS; j++) {
        for (k = 0; k <= j; k++) { // because k is never larger than j
          fprintf(output3, " %.3f", float_edge_stats_bins[i][j][k]);
        }
      }
    }
    fprintf(output3, "\n");
    fclose(output3);
  }
}

// for a given float value, determine which bin it falls into by inspecting the
// upper bounds on each bin
int get_bin(float measure, float *upper_bounds) {
  int i;
  for (i = 0; i < NUM_BINS - 1; i++) {
    if (measure < upper_bounds[i]) {
      return i;
    }
  }
  return NUM_BINS - 1; // i.e. the final bin will contain everything too large
                       // for the others
}
