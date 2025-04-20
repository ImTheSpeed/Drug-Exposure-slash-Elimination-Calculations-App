# Drug Exposre Simul App.

This web application models drug exposure and elimination in the human body using a one-compartment pharmacokinetic model. the app allows users to simulate intravenous (IV) bolus and oral dosing with customizable parameters, including multiple-dose regimens.

## Features
- Test it out on https://math2-drugexpo-proto.streamlit.app/
- Visualize drug concentration over time
- Switch between IV bolus and oral dosing models
- Configure dose, volume of distribution (Vd), and rate constants (k, ka)
- Simulate single or repeated doses with adjustable intervals
- View real-time plots, tabular data, and download results as CSV
- Built-in display of mathematical models and calculated AUC

## Technologies

- Python
- Streamlit
- NumPy, Pandas, Matplotlib

## Getting Started

To run the app locally:

```bash
git clone https://github.com/your-username/pharmacokinetics-app.git
cd pharmacokinetics-app
pip install -r requirements.txt
streamlit run app.py
