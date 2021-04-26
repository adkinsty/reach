import os
import numpy as np
import pandas as pd
from pandas import DataFrame
import os.path
import statsmodels.api as sm
import statsmodels


BEHAV_PATH = os.getcwd()

ALL_SUBS = []  # range(403, 406)

for n in range(300, 600):
    accept_file = BEHAV_PATH + '/' + str(n) + "_rejected.npy"
    reject_file = BEHAV_PATH + '/' + str(n) + "_rejected.npy"
    if os.path.isfile(accept_file) and os.path.isfile(reject_file):
        print(n)
        ALL_SUBS.append(n)
    else:
        continue

print(ALL_SUBS)
ACCEPT = 1
REJECT = 0

POS_OFFER = 0  # The first part of the offer is always positive
NEG_OFFER = 1


def make_loss_frames():
    for sub in ALL_SUBS:

        frame_file = BEHAV_PATH + '/' + str(sub) + "_loss_aversion_frame.csv"
        if os.path.isfile(frame_file):
            print(str(sub) + "exists")
        else:
            decisions = []
            pos = []
            neg = []
            accepted = np.load(BEHAV_PATH + '/' + str(sub) + "_accepted.npy")
            rejected = np.load(BEHAV_PATH + '/' + str(sub) + "_rejected.npy")

            for i, offer in enumerate(accepted):
                decisions.append(ACCEPT)
                pos.append(accepted[i][POS_OFFER])
                neg.append(accepted[i][NEG_OFFER])
            for j, off in enumerate(rejected):
                decisions.append(REJECT)
                pos.append(rejected[j][POS_OFFER])
                neg.append(rejected[j][NEG_OFFER])

            frame_dict = {'decision': decisions, 'gains': pos, 'losses': neg}
            frame = DataFrame(data=frame_dict)
            frame.to_csv(frame_file, sep='\t')


def calc_betas(filename):
    frame = pd.read_table(filename)
    frame['intercept'] = 1.0
    train_cols = frame.columns[2:]
    model = sm.Logit(frame['decision'], frame[train_cols])
    try:
        result = model.fit()
        loss_beta = result.params['losses']
        gain_beta = result.params['gains']
        lambda_ratio = loss_beta / gain_beta
        return loss_beta, gain_beta, lambda_ratio
    except (np.linalg.LinAlgError, statsmodels.tools.sm_exceptions.PerfectSeparationError):
        return np.nan, np.nan, np.nan


def main():
    loss_betas = []
    gain_betas = []
    lambdas = []
    make_loss_frames()
    for sub in ALL_SUBS:
        frame_file = '%s_loss_aversion_frame.csv' % sub
        loss, gain, lam = calc_betas(frame_file)
        print(str(sub) + ":")
        print(str(loss) + " " + str(gain) + " " + str(lam))
        loss_betas.append(loss)
        gain_betas.append(gain)
        lambdas.append(lam)
    l_a_dict = {'subNum': ALL_SUBS, 'gain_betas': gain_betas,
                'loss_betas': loss_betas, 'lambdas': lambdas}
    l_a_frame = DataFrame(data=l_a_dict, index=ALL_SUBS)
    l_a_frame.to_csv('All_lambdas.csv', sep='\t')


if __name__ == '__main__':
    main()
