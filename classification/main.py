from utils.classifier_processor import CentrifugeOutputProcessor, KrakenOutputProcessor, DiamondOutputProcessor, KrakenUniqOutputProcessor

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Process classification output files.")
    parser.add_argument("--input", type=str, help="Path to classification output file.", default=None)
    parser.add_argument("--type", type=str, help="type of classifier used.", default = None, 
                        choices=['centrifuge', 'kraken2', 'diamond', 'kuniq'])
    parser.add_argument("--output_path", type=str, required=True, help="Path to save the final report.")
    parser.add_argument("--nuniq_threshold", type=int, default=1, help="Threshold for unique reads in Centrifuge and Kraken2 output.")

    return parser.parse_args()

def main():
    args = get_args()
    

    if args.type == 'centrifuge':
        centrifuge_processor = CentrifugeOutputProcessor(args.input, nuniq_threshold=args.nuniq_threshold)
        centrifuge_processor.from_file(args.output_path).process().prep_final_report().save(args.output_path)

    if args.type == 'kraken2':
        kraken_processor = KrakenOutputProcessor(args.input, min_uniq_reads=args.nuniq_threshold)
        kraken_processor.from_file(args.output_path).process().prep_final_report().save(args.output_path)

    if args.type == 'diamond':
        diamond_processor = DiamondOutputProcessor(args.input, min_uniq_reads=args.nuniq_threshold)
        diamond_processor.from_file(args.output_path).process().prep_final_report().save(args.output_path)
        
    if args.type == 'kuniq':
        kuniq_processor = KrakenUniqOutputProcessor(args.input, min_uniq_reads=args.nuniq_threshold)
        kuniq_processor.from_file(args.output_path).process().prep_final_report().save(args.output_path)
    
if __name__ == '__main__':
    main()  
    