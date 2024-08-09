import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import pod5
from remora import io, refine_signal_map
import tensorflow as tf
import plotly.graph_objects as go
import plotly

opt_parser = argparse.ArgumentParser()


opt_parser.add_argument(
    "-s",
    "--start_index",
    dest="start_index",
    help="Define a start point of the region of interest on the reference",
    metavar="FILE",
)

opt_parser.add_argument(
    "-e",
    "--end_index",
    dest="end_index",
    help="Define n end point of the region of interest on the reference",
    metavar="FILE",
)

opt_parser.add_argument(
    "-c",
    "--chunk_size",
    dest="chunk_size",
    help="Define the chunk size of the investigatd fragments",
    metavar="FILE",
)

opt_parser.add_argument(
    "-x",
    "--max_seq_length",
    dest="max_seq_length",
    help="Define the maximal sequence length",
    metavar="FILE",
)


opt_parser.add_argument(
    "-r",
    "--reference_path",
    dest="reference_path",
    help="Define a path to the reference",
    metavar="FILE",
)

opt_parser.add_argument(
    "-p",
    "--pod5_path",
    dest="pod5_path",
    help="Define a path to the pod5 of interest",
    metavar="FILE",
)

opt_parser.add_argument(
    "-b",
    "--bam_path",
    dest="bam_path",
    help="Define a path to the bam file of interest",
    metavar="FILE",
)

opt_parser.add_argument(
    "-m",
    "--model",
    dest="model_path",
    help="Define a path to your custom model",
    metavar="FILE",
)

opt_parser.add_argument(
    "-l",
    "--level_table_file",
    dest="level_table_file",
    help="Define a path to ONT's kmer level table. These can be found on github.",
    metavar="FILE",
)


options = opt_parser.parse_args()
start_index = int(options.start_index)
end_index = int(options.end_index)
chunk_size = int(options.chunk_size)
max_seq_length = int(options.max_seq_length)
reference_path = str(options.reference_path)
pod5_path = str(options.pod5_path)
bam_path = str(options.bam_path)
model_path = str(options.model_path)
level_table_file = str(options.level_table_file)


def NN_analyzer(variables, pod5_dr, bam_fh, read_id, sig_map_refiner, model, reference):
    chunck_size = variables[2]
    max_seq_len = variables[3]
    labels = 4
    N_miss = 0
    # ALLERT TO ADJUST!!. maybe there is a problem with the track mode
    reference_track_mod = np.zeros([len(reference), 4])
    start_Index = variables[0]
    for name_id in read_id[variables[0] : variables[1]]:
        pod5_read = pod5_dr.get_read(name_id)
        bam_read = bam_fh.get_first_alignment(name_id)
        start_Index += 1
        print(start_Index)
        seq_resquigle = ""
        position_adjusting = 0
        Error_read = False
        if bam_read.is_reverse:  # correct the signal for forward direction
            flip = False
        else:
            flip = True
        try:
            # /// read data
            read_analysed = io.Read.from_pod5_and_alignment(
                pod5_read, bam_read, reverse_signal=flip
            )
            # // resquigle the data with the refence
            read_analysed.set_refine_signal_mapping(sig_map_refiner, ref_mapping=True)
            start_of_mapping = read_analysed.extract_ref_reg(
                read_analysed.ref_reg.adjust(
                    start_adjust=0, end_adjust=read_analysed.ref_reg.len
                )
            )
            Raw_signal = start_of_mapping.norm_signal
            seq_resquigle = start_of_mapping.seq
            start_end_resquigle = start_of_mapping.seq_to_sig_map
            # /// check if the modification position has to be adjusted ///
            position_adjusting = start_of_mapping.ref_reg.start
        except:
            print("error")
            position_adjusting = 0
            seq_resquigle = ""
            Error_read = True
        if not Error_read:
            base_dict = {"A": 1, "C": 2, "G": 3, "T": 4}
            bases_onehot = np.zeros([len(Raw_signal), 4 + 1])
            for k in range(len(seq_resquigle)):
                start_resq = start_end_resquigle[k]
                bases_onehot[start_resq, base_dict[seq_resquigle[k]]] = 1
            try:
                N_segments = int(len(Raw_signal) / chunck_size)
                Input_1 = np.zeros(
                    [N_segments + 1, chunck_size]
                )  # initialize the first input of the NN
                Input_2 = np.zeros(
                    [N_segments + 1, max_seq_len, labels]
                )  # initialize the second input of the NN
                for k in range(N_segments):
                    start = k * chunck_size
                    Input_1[k] = Raw_signal[start : start + chunck_size]
                    window_onehot = bases_onehot[start : start + chunck_size, :]
                    probe = np.argmax(window_onehot, axis=-1)
                    probe = probe[probe != 0]
                    probe = probe - 1
                    for kk in range(len(probe)):
                        Input_2[k, kk, probe[kk]] = 1
                # find the number of point not overlapping
                not_overlaping_last_seg = len(Raw_signal) - (start + chunck_size)
                # the extention to +1 is for keeping the full dimention of the output
                Input_1[N_segments] = Raw_signal[-chunck_size:]
                Additional_window = bases_onehot[-chunck_size:, :]
                probe = np.argmax(Additional_window, axis=-1)
                probe = probe[probe != 0]
                probe = probe - 1
                for kk in range(len(probe)):
                    Input_2[N_segments, kk, probe[kk]] = 1
                # probe the overlapping bases for the last segment
                Window_overlap = bases_onehot[-chunck_size:-not_overlaping_last_seg, :]
                seq_overlap = np.zeros([Window_overlap.shape[0], 4])
                probe = np.argmax(Window_overlap, axis=-1)
                probe = probe[probe != 0]
                probe = probe - 1
                for kk in range(len(probe)):
                    seq_overlap[kk, probe[kk]] = 1
                seq_overlap = np.sum(seq_overlap, axis=1)
                seq_overlap = np.where(seq_overlap > 0.5)[0]
                len_overlap = len(seq_overlap)
                Input_1 = np.expand_dims(Input_1, axis=-1)
                Input_2 = np.expand_dims(Input_2, axis=-1)
                X_total = {"Input_1": Input_1, "Input_2": Input_2}
                # analyze the read with the NN
                prediction = model.predict(X_total)
                # reconstruct the final output removing the null part of the predictions
                Final_seq_binary = []
                for kk in range(N_segments):  #
                    full_position = np.sum(prediction[kk], axis=1)
                    full_position = np.where(full_position > 0.5)[0]
                    real_part = np.argmax(prediction[kk, : len(full_position)], axis=-1)
                    Final_seq_binary = np.concatenate(
                        (Final_seq_binary, real_part), axis=0
                    )
                full_position = np.sum(prediction[N_segments], axis=1)
                full_position = np.where(full_position > 0.5)[0]
                real_part = np.argmax(
                    prediction[N_segments, : len(full_position)], axis=-1
                )
                not_overlaping_part = real_part[len_overlap:]
                Final_seq_binary = np.concatenate(
                    (Final_seq_binary, not_overlaping_part), axis=0
                )
                if (len(Final_seq_binary) - len(seq_resquigle)) != 0:
                    print("ops, missmatch:", len(Final_seq_binary) - len(seq_resquigle))
                    N_miss += 1
                else:
                    where_mod = np.where(Final_seq_binary >= 1)[0]
                    modific_detec = np.zeros(len(where_mod))
                    for j in range(len(where_mod)):
                        modific_detec[j] = Final_seq_binary[where_mod[j]]

                    if len(modific_detec) > 1:
                        for n in range(len(modific_detec)):
                            mod_probe_position = where_mod[n]
                            mod_probe_predicted = modific_detec[n]
                            reference_track_mod[
                                int(mod_probe_position) + int(position_adjusting),
                                int(mod_probe_predicted - 1),
                            ] += 1
                    else:
                        mod_probe_position = where_mod[0]
                        mod_probe_predicted = modific_detec[0]
                        reference_track_mod[
                            int(mod_probe_position) + int(position_adjusting),
                            int(mod_probe_predicted - 1),
                        ] += 1
            except:
                print("No mod detected")
                N_miss += 1

    return (reference_track_mod) / (np.abs(variables[1] - variables[0]))


def Analysis_Neural_network(
    start_index: int,
    end_index: int,
    chunk_size: int,
    max_seq_length: int,
    reference_path: str,
    pod5_path: str,
    bam_path: str,
    model_path: str,
    level_table_file: str,
):

    vars_entries = [start_index, end_index, chunk_size, max_seq_length]
    pod5_dr = pod5.DatasetReader(pod5_path)
    bam_fh = io.ReadIndexedBam(bam_path)
    read_id = bam_fh.read_ids
    NN_model = tf.keras.models.load_model(model_path)
    sig_map_refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=level_table_file,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )
    Variables = (vars_entries[0], vars_entries[1], vars_entries[2], vars_entries[3])

    reference = open(reference_path)
    reference = reference.read()

    Analysis_NN = NN_analyzer(
        Variables, pod5_dr, bam_fh, read_id, sig_map_refiner, NN_model, reference
    )

    x_axis = np.arange(0, Analysis_NN.shape[0], 1)
    print(Analysis_NN)

    layout = go.Layout(height=800)
    fig = go.Figure(layout=layout)

    hover_text = [f"Gm: {y_val}" for y_val in Analysis_NN[:, 0]]
    fig.add_trace(
        go.Scatter(
            x=x_axis,
            y=Analysis_NN[:, 0],
            mode="lines+markers",
            line=dict(color="rgba(72,99,156,1)"),
            hovertext=hover_text,
            hoverinfo = "text"
        )
    )

    hover_text = [f"$m^6A$: {y_val}" for y_val in Analysis_NN[:, 1]]
    fig.add_trace(
        go.Scatter(
            x=x_axis,
            y=Analysis_NN[:, 1],
            mode="lines+markers",
            line=dict(color="rgba(95,62,139,1)"),
            hovertext=hover_text,
            hoverinfo = "text"
        )
    )

    hover_text = [f"Ino: {y_val}" for y_val in Analysis_NN[:, 2]]
    fig.add_trace(
        go.Scatter(
            x=x_axis,
            y=Analysis_NN[:, 2],
            mode="lines+markers",            
            line=dict(color="rgba(132,73,133,1)"),
            hovertext=hover_text,
            hoverinfo = "text"
        )
    )

    hover_text = [f"Psi: {y_val}" for y_val in Analysis_NN[:, 3]]
    fig.add_trace(
        go.Scatter(
            x=x_axis,
            y=Analysis_NN[:, 3],
            mode="lines+markers",            
            line=dict(color="rgba(247,153,110,1)"),
            hovertext=hover_text,
            hoverinfo = "text"
        )
    )

    fig.update_layout(
        xaxis=dict(title="Position on reference", gridcolor="white", tickformat="d"),
        yaxis=dict(
            title="Modification frequency",
            gridcolor="white",
            zeroline=True,
            zerolinecolor="black",
            zerolinewidth=0.1,
            tickformat="d",
        ),
        plot_bgcolor="rgba(0,0,0,0)",
    )

    plotly.io.write_html(fig, "./detected_modifications.html")


Analysis_Neural_network(
    start_index,
    end_index,
    chunk_size,
    max_seq_length,
    reference_path,
    pod5_path,
    bam_path,
    model_path,
    level_table_file,
)
