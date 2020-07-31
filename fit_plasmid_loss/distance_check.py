import pandas as pd

def drop_repeats():
    df = pd.read_csv('./GS_plasmid_loss.csv')
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    df.drop_duplicates(subset=['plasmid_name', 'w_val', 'init_Q_N', 'PCN'], keep=False, inplace=True)
    
    df.to_csv('./GS_plasmid_loss_2.csv')

def main():
    df = pd.read_csv('./GS_plasmid_loss.csv')
    plasmid_names = df['plasmid_name'].unique()
    w_vals = df['w_val'].unique()

    total_min_distance =1e6
    best_w_val = None
    
    for w in w_vals:
        total_distance = 0
        w_sub_df = df.loc[df['w_val'] == w ]
        for plasmid in plasmid_names:
            p_sub_df = w_sub_df.loc[w_sub_df['plasmid_name'] == plasmid].reset_index()
            
            if w == 30.0:
                print(p_sub_df.iloc[p_sub_df['distance'].idxmin()])
            
            min_dist = p_sub_df.iloc[p_sub_df['distance'].idxmin()]['distance']
            total_distance += min_dist

        if total_distance < total_min_distance:
            total_min_distance = total_min_distance
            best_w_val = w

        print(w, total_distance)
        print("")

    print("")
    print(best_w_val)

if __name__== "__main__":
    drop_repeats()