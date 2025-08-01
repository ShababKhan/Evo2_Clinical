"""
Treatment Outcome Visualization Module.

This module provides visualization components for treatment outcomes,
including response predictions, survival curves, and adverse event risk.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Union, Tuple, Any
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import io
import base64
from matplotlib.colors import LinearSegmentedColormap


class TreatmentOutcomeVisualizer:
    """Class for visualizing treatment outcomes and predictions."""
    
    def __init__(self, theme: str = "light"):
        """
        Initialize the treatment outcome visualizer.
        
        Args:
            theme (str): Color theme for visualizations ("light" or "dark")
        """
        self.theme = theme
        self._setup_styles()
    
    def _setup_styles(self):
        """Set up visualization styles based on theme."""
        if self.theme == "dark":
            self.colors = {
                "background": "#1E1E1E",
                "text": "#FFFFFF",
                "treatment_types": {
                    "radiation": "#BB86FC",
                    "chemotherapy": "#03DAC6",
                    "immunotherapy": "#CF6679",
                    "targeted": "#FFB74D",
                    "other": "#BBBBBB"
                },
                "outcomes": {
                    "response": "#64FFDA", 
                    "survival": "#4FC3F7",
                    "adverse": "#EF5350"
                }
            }
            plt.style.use("dark_background")
        else:
            self.colors = {
                "background": "#FFFFFF",
                "text": "#000000",
                "treatment_types": {
                    "radiation": "#7B1FA2",
                    "chemotherapy": "#00796B",
                    "immunotherapy": "#C62828",
                    "targeted": "#EF6C00",
                    "other": "#757575"
                },
                "outcomes": {
                    "response": "#00695C", 
                    "survival": "#0288D1",
                    "adverse": "#D32F2F"
                }
            }
            plt.style.use("default")
    
    def plot_treatment_comparison(
        self,
        treatments_data: Union[pd.DataFrame, Dict[str, Dict[str, float]]],
        response_col: str = "response_probability",
        survival_col: str = "survival_probability",
        adverse_col: str = "adverse_event_risk",
        treatment_col: str = "treatment",
        title: str = "Treatment Outcome Comparison",
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot comparison of treatment outcomes.
        
        Args:
            treatments_data: DataFrame with treatment data or dictionary mapping treatments to outcome metrics
            response_col (str): Column with response probabilities
            survival_col (str): Column with survival probabilities
            adverse_col (str): Column with adverse event risks
            treatment_col (str): Column with treatment names
            title (str): Plot title
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        # Convert dictionary to DataFrame if necessary
        if isinstance(treatments_data, dict):
            data_list = []
            for treatment, outcomes in treatments_data.items():
                row = {"treatment": treatment}
                row.update(outcomes)
                data_list.append(row)
            df = pd.DataFrame(data_list)
        else:
            df = treatments_data
        
        # Check for required columns
        required_cols = [treatment_col, response_col, survival_col, adverse_col]
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"Warning: Missing columns: {missing_cols}")
            return None
        
        # For interactive plots with Plotly
        if interactive:
            # Create a subplot with shared x-axis
            fig = make_subplots(
                rows=3, 
                cols=1,
                subplot_titles=("Response Probability", "Survival Probability", "Adverse Event Risk"),
                shared_xaxes=True,
                vertical_spacing=0.1
            )
            
            # Sort treatments by response probability for consistent ordering
            sorted_df = df.sort_values(response_col, ascending=False).reset_index(drop=True)
            treatments = sorted_df[treatment_col]
            
            # Add traces for each outcome
            # Response probability
            fig.add_trace(
                go.Bar(
                    x=treatments,
                    y=sorted_df[response_col],
                    marker_color=self.colors["outcomes"]["response"],
                    name="Response",
                    hovertemplate="%{x}: %{y:.3f}"
                ),
                row=1, col=1
            )
            
            # Survival probability
            fig.add_trace(
                go.Bar(
                    x=treatments,
                    y=sorted_df[survival_col],
                    marker_color=self.colors["outcomes"]["survival"],
                    name="Survival",
                    hovertemplate="%{x}: %{y:.3f}"
                ),
                row=2, col=1
            )
            
            # Adverse event risk
            fig.add_trace(
                go.Bar(
                    x=treatments,
                    y=sorted_df[adverse_col],
                    marker_color=self.colors["outcomes"]["adverse"],
                    name="Adverse Event",
                    hovertemplate="%{x}: %{y:.3f}"
                ),
                row=3, col=1
            )
            
            # Update layout
            fig.update_layout(
                title=title,
                height=600,
                showlegend=False,
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            
            # Update y-axis ranges
            fig.update_yaxes(range=[0, 1], title_text="Probability", row=1, col=1)
            fig.update_yaxes(range=[0, 1], title_text="Probability", row=2, col=1)
            fig.update_yaxes(range=[0, 1], title_text="Risk", row=3, col=1)
            
            return fig
        else:
            # For static plots with Matplotlib
            fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
            
            # Sort treatments by response probability for consistent ordering
            sorted_df = df.sort_values(response_col, ascending=False).reset_index(drop=True)
            treatments = sorted_df[treatment_col]
            
            # Response probability
            axs[0].bar(
                treatments, 
                sorted_df[response_col], 
                color=self.colors["outcomes"]["response"]
            )
            axs[0].set_title("Response Probability")
            axs[0].set_ylabel("Probability")
            axs[0].set_ylim(0, 1)
            
            # Survival probability
            axs[1].bar(
                treatments, 
                sorted_df[survival_col], 
                color=self.colors["outcomes"]["survival"]
            )
            axs[1].set_title("Survival Probability")
            axs[1].set_ylabel("Probability")
            axs[1].set_ylim(0, 1)
            
            # Adverse event risk
            axs[2].bar(
                treatments, 
                sorted_df[adverse_col], 
                color=self.colors["outcomes"]["adverse"]
            )
            axs[2].set_title("Adverse Event Risk")
            axs[2].set_xlabel("Treatment")
            axs[2].set_ylabel("Risk")
            axs[2].set_ylim(0, 1)
            
            # Rotate treatment labels
            plt.setp(axs[2].xaxis.get_majorticklabels(), rotation=45, ha="right")
            
            plt.suptitle(title, fontsize=16)
            plt.tight_layout()
            
            return fig
    
    def plot_detailed_adverse_events(
        self,
        adverse_events_data: Union[pd.DataFrame, Dict[str, Dict[str, float]]],
        treatment_col: str = "treatment",
        overall_col: str = "overall",
        title: str = "Adverse Event Risk by Treatment Type",
        max_events: int = 10,
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot detailed adverse event risks by treatment type.
        
        Args:
            adverse_events_data: DataFrame with adverse event data or dictionary 
                               mapping treatments to event risks
            treatment_col (str): Column with treatment names
            overall_col (str): Column/key for overall adverse event risk
            title (str): Plot title
            max_events (int): Maximum number of adverse events to show
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        # Convert dictionary to DataFrame if necessary
        if isinstance(adverse_events_data, dict):
            data_list = []
            for treatment, events in adverse_events_data.items():
                for event, risk in events.items():
                    if event != overall_col:  # Skip overall risk for this plot
                        data_list.append({
                            "treatment": treatment,
                            "event": event,
                            "risk": risk
                        })
            df = pd.DataFrame(data_list)
        else:
            # Reshape DataFrame to long format
            df_long = pd.melt(
                adverse_events_data,
                id_vars=[treatment_col],
                var_name="event",
                value_name="risk"
            )
            df = df_long[df_long["event"] != overall_col]  # Skip overall risk
        
        # Ensure we have the required columns
        if not set(["treatment", "event", "risk"]).issubset(df.columns):
            print("Warning: Required columns missing")
            return None
        
        # Limit to max_events most common adverse events
        top_events = df.groupby("event")["risk"].mean().sort_values(ascending=False).head(max_events).index
        df = df[df["event"].isin(top_events)]
        
        # For interactive plots
        if interactive:
            # Create a heatmap with Plotly
            pivot_df = df.pivot_table(values="risk", index="event", columns="treatment", aggfunc="max")
            
            # Map treatment types to colors if available
            colorscale = [
                [0, "green"],     # Low risk
                [0.5, "yellow"],  # Medium risk
                [1, "red"]        # High risk
            ]
            
            fig = px.imshow(
                pivot_df,
                color_continuous_scale=colorscale,
                labels=dict(x="Treatment", y="Adverse Event", color="Risk"),
                title=title
            )
            
            # Add text annotations
            for i, event in enumerate(pivot_df.index):
                for j, treatment in enumerate(pivot_df.columns):
                    risk = pivot_df.iloc[i, j]
                    fig.add_annotation(
                        x=treatment,
                        y=event,
                        text=f"{risk:.2f}",
                        showarrow=False,
                        font=dict(
                            color="black" if risk < 0.5 else "white"
                        )
                    )
            
            fig.update_layout(
                width=800,
                height=600,
                template="plotly_white" if self.theme == "light" else "plotly_dark"
            )
            
            return fig
        else:
            # For static plots with Matplotlib
            plt.figure(figsize=(12, 8))
            
            # Create pivot table
            pivot_df = df.pivot_table(values="risk", index="event", columns="treatment", aggfunc="max")
            
            # Create a custom colormap: green to yellow to red
            cmap = LinearSegmentedColormap.from_list("risk_cmap", ["green", "yellow", "red"])
            
            # Plot heatmap
            ax = sns.heatmap(
                pivot_df,
                cmap=cmap,
                annot=True,
                fmt=".2f",
                linewidths=0.5
            )
            
            plt.title(title)
            plt.ylabel("Adverse Event")
            plt.xlabel("Treatment")
            plt.tight_layout()
            
            return plt.gcf()
    
    def plot_survival_curve(
        self,
        time_points: Union[List[float], np.ndarray],
        survival_probs: Union[Dict[str, List[float]], Dict[str, np.ndarray]],
        confidence_intervals: Optional[Dict[str, List[Tuple[float, float]]]] = None,
        title: str = "Predicted Survival Curve",
        xlabel: str = "Time (months)",
        ylabel: str = "Survival Probability",
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot Kaplan-Meier style survival curves.
        
        Args:
            time_points (List[float]): Time points for x-axis
            survival_probs (Dict[str, List[float]]): Dictionary mapping 
                                                     treatment names to survival probabilities
            confidence_intervals (Dict[str, List[Tuple[float, float]]], optional): 
                                                Dictionary mapping treatment names to 
                                                confidence intervals for each time point
            title (str): Plot title
            xlabel (str): X-axis label
            ylabel (str): Y-axis label
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if not survival_probs:
            return None
        
        # For interactive plots with Plotly
        if interactive:
            fig = go.Figure()
            
            # Add a survival curve for each treatment
            for i, (treatment, probs) in enumerate(survival_probs.items()):
                # Get a color for this treatment type
                if treatment.lower() in self.colors["treatment_types"]:
                    color = self.colors["treatment_types"][treatment.lower()]
                else:
                    # Use the 'other' color or cycle through a default set
                    color = self.colors["treatment_types"].get("other", "#757575")
                
                # Add the main survival curve
                fig.add_trace(go.Scatter(
                    x=time_points,
                    y=probs,
                    mode='lines',
                    name=treatment,
                    line=dict(color=color)
                ))
                
                # Add confidence intervals if available
                if confidence_intervals and treatment in confidence_intervals:
                    ci_data = confidence_intervals[treatment]
                    lower_bound = [ci[0] for ci in ci_data]
                    upper_bound = [ci[1] for ci in ci_data]
                    
                    # Add the confidence interval as a filled area
                    fig.add_trace(go.Scatter(
                        x=list(time_points) + list(reversed(time_points)),
                        y=list(upper_bound) + list(reversed(lower_bound)),
                        fill='toself',
                        fillcolor=f'rgba({int(color[1:3], 16)}, {int(color[3:5], 16)}, {int(color[5:7], 16)}, 0.2)',
                        line=dict(color='rgba(0,0,0,0)'),
                        name=f'{treatment} CI',
                        showlegend=False
                    ))
            
            # Update layout
            fig.update_layout(
                title=title,
                xaxis_title=xlabel,
                yaxis_title=ylabel,
                yaxis=dict(range=[0, 1.05]),
                template="plotly_white" if self.theme == "light" else "plotly_dark",
                hovermode="x unified"
            )
            
            return fig
        else:
            # For static plots with Matplotlib
            plt.figure(figsize=(10, 6))
            
            # Add a survival curve for each treatment
            for treatment, probs in survival_probs.items():
                # Get a color for this treatment type
                if treatment.lower() in self.colors["treatment_types"]:
                    color = self.colors["treatment_types"][treatment.lower()]
                else:
                    color = self.colors["treatment_types"].get("other", "#757575")
                
                # Plot the main survival curve
                plt.step(time_points, probs, where='post', label=treatment, color=color)
                
                # Add confidence intervals if available
                if confidence_intervals and treatment in confidence_intervals:
                    ci_data = confidence_intervals[treatment]
                    lower_bound = [ci[0] for ci in ci_data]
                    upper_bound = [ci[1] for ci in ci_data]
                    
                    plt.fill_between(
                        time_points,
                        lower_bound,
                        upper_bound,
                        alpha=0.2,
                        step='post',
                        color=color
                    )
            
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.ylim(0, 1.05)
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.legend()
            plt.tight_layout()
            
            return plt.gcf()
    
    def plot_dose_response(
        self,
        doses: Union[List[float], np.ndarray],
        responses: Union[Dict[str, List[float]], Dict[str, np.ndarray]],
        title: str = "Dose-Response Curve",
        xlabel: str = "Dose",
        ylabel: str = "Response",
        show_ec50: bool = True,
        interactive: bool = True
    ) -> Union[go.Figure, plt.Figure]:
        """
        Plot dose-response curves for different treatments or variants.
        
        Args:
            doses (List[float]): Dose levels
            responses (Dict[str, List[float]]): Dictionary mapping names to response values
            title (str): Plot title
            xlabel (str): X-axis label
            ylabel (str): Y-axis label
            show_ec50 (bool): Whether to highlight EC50 values
            interactive (bool): Whether to use interactive Plotly or static Matplotlib
            
        Returns:
            Union[go.Figure, plt.Figure]: Plotly or Matplotlib figure
        """
        if not responses:
            return None
        
        # Calculate EC50 values if requested
        ec50_values = {}
        if show_ec50:
            for name, resp in responses.items():
                # Find the dose closest to 50% response
                if max(resp) > 0.5 > min(resp):
                    # Simple linear interpolation
                    for i in range(len(doses) - 1):
                        if resp[i] <= 0.5 <= resp[i+1] or resp[i] >= 0.5 >= resp[i+1]:
                            # Linear interpolation to find EC50
                            if resp[i] != resp[i+1]:  # Avoid division by zero
                                ec50 = doses[i] + (doses[i+1] - doses[i]) * (0.5 - resp[i]) / (resp[i+1] - resp[i])
                                ec50_values[name] = ec50
                            break
        
        # For interactive plots with Plotly
        if interactive:
            fig = go.Figure()
            
            # Add a dose-response curve for each item
            for i, (name, resp) in enumerate(responses.items()):
                # Get a color for this curve
                color_key = name.lower() if name.lower() in self.colors["treatment_types"] else "other"
                color = self.colors["treatment_types"].get(color_key, "#757575")
                
                # Add the dose-response curve
                fig.add_trace(go.Scatter(
                    x=doses,
                    y=resp,
                    mode='lines+markers',
                    name=name,
                    line=dict(color=color)
                ))
                
                # Add EC50 marker if available
                if name in ec50_values:
                    ec50 = ec50_values[name]
                    fig.add_trace(go.Scatter(
                        x=[ec50],
                        y=[0.5],
                        mode='markers',
                        marker=dict(
                            size=10,
                            color=color,
                            symbol='diamond',
                            line=dict(
                                color='black',
                                width=2
                            )
                        ),
                        name=f"{name} EC50",
                        text=f"EC50: {ec50:.3f}",
                        hoverinfo='text'
                    ))
                    
                    # Add EC50 reference lines
                    fig.add_shape(
                        type="line",
                        x0=ec50,
                        y0=0,
                        x1=ec50,
                        y1=0.5,
                        line=dict(color=color, dash="dash", width=1),
                        opacity=0.7
                    )
                    fig.add_shape(
                        type="line",
                        x0=min(doses),
                        y0=0.5,
                        x1=ec50,
                        y1=0.5,
                        line=dict(color=color, dash="dash", width=1),
                        opacity=0.7
                    )
            
            # Update layout
            fig.update_layout(
                title=title,
                xaxis_title=xlabel,
                yaxis_title=ylabel,
                template="plotly_white" if self.theme == "light" else "plotly_dark",
                xaxis_type="log" if (max(doses) / min(doses) > 100) else "linear"
            )
            
            return fig
        else:
            # For static plots with Matplotlib
            plt.figure(figsize=(10, 6))
            
            # Add a dose-response curve for each item
            for name, resp in responses.items():
                # Get a color for this curve
                color_key = name.lower() if name.lower() in self.colors["treatment_types"] else "other"
                color = self.colors["treatment_types"].get(color_key, "#757575")
                
                # Plot the dose-response curve
                plt.plot(doses, resp, 'o-', label=name, color=color)
                
                # Add EC50 marker if available
                if name in ec50_values:
                    ec50 = ec50_values[name]
                    plt.scatter([ec50], [0.5], color=color, marker='D', s=100, edgecolors='black')
                    
                    # Add EC50 reference lines
                    plt.axvline(x=ec50, color=color, linestyle='--', alpha=0.5)
                    plt.axhline(y=0.5, color=color, linestyle='--', alpha=0.5)
                    
                    # Add EC50 label
                    plt.text(
                        ec50, 
                        0.55, 
                        f"EC50: {ec50:.3f}", 
                        color=color, 
                        fontweight='bold',
                        horizontalalignment='center'
                    )
            
            # Use log scale if dose range is large
            if max(doses) / min(doses) > 100:
                plt.xscale('log')
            
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.legend()
            plt.tight_layout()
            
            return plt.gcf()
    
    def fig_to_base64(self, fig: Union[go.Figure, plt.Figure]) -> str:
        """
        Convert a figure to base64 string for embedding in HTML/Markdown.
        
        Args:
            fig: Matplotlib or Plotly figure object
            
        Returns:
            str: Base64 encoded image string
        """
        if fig is None:
            return ""
            
        buffer = io.BytesIO()
        
        # Handle different figure types
        if isinstance(fig, go.Figure):
            fig.write_image(buffer, format="png")
        else:
            # Matplotlib figure
            fig.savefig(buffer, format="png", bbox_inches="tight")
            plt.close(fig)
        
        # Get the image as base64 string
        buffer.seek(0)
        image_base64 = base64.b64encode(buffer.read()).decode("utf-8")
        return f"data:image/png;base64,{image_base64}"