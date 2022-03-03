import os
import argparse

def verification_bsub_output(file, rerun_command='', add_comment=''):
    print('QUALITY ANALYSIS')
    with open(file, 'r') as f:
        out = f.readlines()
    if 'Successfully completed' in out:
        print(f'\n Successfully completed - {add_comment}')
    elif 'TERM_MEM' in out:
        if rerun_command:
            os.system(rerun_command)
        print(f'\n Rerunning - {add_comment} - {rerun_command}')
    else:
        print(f'\n Unknown error - {add_comment}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verification bsub outputs')
    parser.add_argument('--file', type=str, help='output file')
    parser.add_argument('-r', '--rerun_command', type=str, help='Rerun command in case of LSF shortage')
    parser.add_argument('-a', '--add_comment', type=str, help='Additional comment')
    args = parser.parse_args()
    verification_bsub_output(**vars(args))
