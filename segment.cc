
#include "segment.h"

using namespace std;

void segmentChannels(ATOM_NETWORK *uc, vector<FEATURE> *channel_list,
                     ostream &output) {
  int i;
  vector<double> segment_list;
  if (channel_list->size() > 0) {
    for (i = 0; i < (int)channel_list->size(); i++) {
      channel_list->at(i).segmentChannel(uc);
      segment_list.insert(segment_list.begin(),
                          channel_list->at(i).segment_max_r.begin(),
                          channel_list->at(i).segment_max_r.end());
    }
  }

  output << segment_list.size() << "\n";
  for (i = 0; i < (int)segment_list.size(); i++) {
    output << segment_list[i] << "\n";
  }
}
/* end of segmentChannels */

/* This function runs segmentation of channels and provides two vectors: seg_ids
 * - with segment IDs for each node in voronoi network and segmaxr with max
 * radius assigned to each segment (if featflag set to 0) */
/* the same function is used to do this for features (featflag=1) */
void segmentChannels_forHolograms(ATOM_NETWORK *uc,
                                  vector<FEATURE> *channel_list,
                                  vector<int> *seg_ids, vector<double> *segmaxr,
                                  int featflag) {
  int i;
  int current_segment = 1;
  vector<double> segment_list;

  if (channel_list->size() > 0) {
    for (i = 0; i < (int)channel_list->size(); i++) {
      channel_list->at(i).segmentChannel(uc);
      channel_list->at(i).update_node_segmentinfo(seg_ids, &current_segment,
                                                  featflag);

      if (featflag == 0)
        segmaxr->insert(segmaxr->end(),
                        channel_list->at(i).segment_max_r.begin(),
                        channel_list->at(i).segment_max_r.end());
      else
        segmaxr->insert(segmaxr->end(),
                        channel_list->at(i).feature_max_r.begin(),
                        channel_list->at(i).feature_max_r.end());

      //      segment_list.insert(segment_list.begin(),channel_list->at(i).segment_max_r.begin(),channel_list->at(i).segment_max_r.end());
    }
  }
}

/* end of segmentChannels */

/* This function save active segment files to run simluations for particular
 * segments */

void segmentChannels_saveSegments(ATOM_NETWORK *uc, VORONOI_NETWORK *vornet,
                                  vector<int> *accessInfo,
                                  vector<double> *segment_radii_vector,
                                  int num_segments, char *name,
                                  const char *feattype) {
  int i, j, no;

  char filename_out[256];
  fstream file_out;
  // double x,y,z;
  Point xyzpt, abcpt;

  for (i = 0; i < num_segments; i++) {
    //   sprintf(filename_out,"%s-segment_%d.active",name,(i+1));
    sprintf(filename_out, "%s-%s_%d.active", name, feattype, (i + 1));

    file_out.open(filename_out, fstream::out);
    no = 0;
    for (j = 0; j < (int)(accessInfo->size()); j++) {
      if (accessInfo->at(j) == (i + 1))
        no++;
    };

    file_out << no << "\n";

    for (j = 0; j < (int)(accessInfo->size()); j++) {
      if (accessInfo->at(j) == (i + 1)) {
        xyzpt[0] = vornet->nodes.at(j).x;
        xyzpt[1] = vornet->nodes.at(j).y;
        xyzpt[2] = vornet->nodes.at(j).z;
        abcpt = uc->xyz_to_abc(xyzpt);
        file_out << abcpt[0] << "   " << abcpt[1] << "   " << abcpt[2]
                 << "    ";
        file_out << vornet->nodes.at(j).rad_stat_sphere * SCALE_SPHERE << "\n";
      }
    }
    file_out.flush();
    file_out.close();
  }
}
/* end of saveSegments */
