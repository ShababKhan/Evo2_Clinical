"""
Machine Learning Pipeline for Evo2_Clinical.

This module provides machine learning capabilities for predicting outcomes based on
variant data, including prediction of treatment responses and radiation/drug-induced
lung injury.
"""

import os
import pickle
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union, Any, Tuple
import matplotlib.pyplot as plt
import seaborn as sns

# Machine learning imports
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, roc_curve, confusion_matrix, classification_report
)
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.decomposition import PCA
import shap

# Import project config
from evo2_clinical.config import config


class FeatureExtractor:
    """Extract features from variant and lncRNA data for machine learning."""
    
    def __init__(self, config_obj=None):
        """
        Initialize the feature extractor.
        
        Args:
            config_obj: Configuration object (default: use global config)
        """
        self.config = config_obj if config_obj is not None else config
        self.logger = logging.getLogger(__name__)
    
    def extract_variant_features(
        self,
        variants_df: pd.DataFrame,
        functional_scores: Optional[Dict[str, float]] = None
    ) -> pd.DataFrame:
        """
        Extract features from variant data.
        
        Args:
            variants_df: DataFrame containing variant data
            functional_scores: Dictionary mapping variant IDs to functional scores
            
        Returns:
            DataFrame with extracted features
        """
        if variants_df.empty:
            self.logger.warning("Empty variants DataFrame provided")
            return pd.DataFrame()
        
        # Create a copy to avoid modifying the original
        features_df = variants_df.copy()
        
        # Extract basic variant features
        if 'REF' in features_df.columns and 'ALT' in features_df.columns:
            # Calculate variant length difference
            features_df['length_diff'] = features_df['ALT'].str.len() - features_df['REF'].str.len()
            
            # Determine variant type
            features_df['is_snp'] = (features_df['REF'].str.len() == 1) & (features_df['ALT'].str.len() == 1)
            features_df['is_insertion'] = (features_df['REF'].str.len() < features_df['ALT'].str.len())
            features_df['is_deletion'] = (features_df['REF'].str.len() > features_df['ALT'].str.len())
        
        # Incorporate functional scores if provided
        if functional_scores is not None:
            # Create a variant ID column if not present
            if 'variant_id' not in features_df.columns and all(col in features_df.columns for col in ['CHROM', 'POS', 'REF', 'ALT']):
                features_df['variant_id'] = features_df['CHROM'].astype(str) + '_' + \
                                            features_df['POS'].astype(str) + '_' + \
                                            features_df['REF'] + '_' + \
                                            features_df['ALT']
            
            # Add functional scores
            if 'variant_id' in features_df.columns:
                features_df['functional_score'] = features_df['variant_id'].map(functional_scores)
        
        # Process any gene annotations
        if 'GENE' in features_df.columns:
            # Placeholder for gene-specific features
            pass
        
        # Process any consequences
        if 'CONSEQUENCE' in features_df.columns:
            # One-hot encode consequences
            consequences = pd.get_dummies(features_df['CONSEQUENCE'], prefix='consequence')
            features_df = pd.concat([features_df, consequences], axis=1)
        
        # Drop non-numeric and identifier columns for ML
        non_feature_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 
                           'variant_id', 'CONSEQUENCE', 'GENE']
        feature_cols = [col for col in features_df.columns if col not in non_feature_cols]
        
        # Return only feature columns
        return features_df[feature_cols].fillna(0)
    
    def extract_lncrna_features(
        self,
        lncrna_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Extract features from lncRNA data.
        
        Args:
            lncrna_df: DataFrame containing lncRNA data
            
        Returns:
            DataFrame with extracted features
        """
        if lncrna_df.empty:
            self.logger.warning("Empty lncRNA DataFrame provided")
            return pd.DataFrame()
        
        # Create a copy to avoid modifying the original
        features_df = lncrna_df.copy()
        
        # Extract basic lncRNA features
        # Placeholder for lncRNA-specific feature extraction
        
        # Drop non-numeric and identifier columns for ML
        non_feature_cols = ['lncrna_id', 'name', 'gene_id', 'transcript_id']
        feature_cols = [col for col in features_df.columns if col not in non_feature_cols]
        
        # Return only feature columns
        return features_df[feature_cols].fillna(0)
    
    def combine_features(
        self,
        variant_features: pd.DataFrame,
        lncrna_features: Optional[pd.DataFrame] = None,
        clinical_data: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Combine different feature sets into a single feature matrix.
        
        Args:
            variant_features: DataFrame with variant features
            lncrna_features: DataFrame with lncRNA features (optional)
            clinical_data: DataFrame with clinical data (optional)
            
        Returns:
            DataFrame with combined features
        """
        feature_dfs = [variant_features]
        
        if lncrna_features is not None and not lncrna_features.empty:
            # Ensure there's a common index or join key
            # For now, we'll assume the index can be used for joining
            feature_dfs.append(lncrna_features)
        
        if clinical_data is not None and not clinical_data.empty:
            # Process clinical data to extract relevant features
            # For now, we'll assume the clinical data is already processed
            feature_dfs.append(clinical_data)
        
        # Combine all feature DataFrames
        combined_features = pd.concat(feature_dfs, axis=1)
        
        # Handle any duplicate columns
        combined_features = combined_features.loc[:, ~combined_features.columns.duplicated()]
        
        return combined_features


class ModelTrainer:
    """Train and evaluate machine learning models for outcome prediction."""
    
    def __init__(self, config_obj=None):
        """
        Initialize the model trainer.
        
        Args:
            config_obj: Configuration object (default: use global config)
        """
        self.config = config_obj if config_obj is not None else config
        self.logger = logging.getLogger(__name__)
        self.models = {}
        self.model_metrics = {}
        self.feature_importances = {}
        self.shap_values = {}
    
    def prepare_data(
        self,
        features: pd.DataFrame,
        target: pd.Series,
        test_size: float = 0.2,
        random_state: int = 42
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series]:
        """
        Prepare data for model training by splitting into train and test sets.
        
        Args:
            features: Feature DataFrame
            target: Target Series
            test_size: Proportion of data to use for testing
            random_state: Random seed for reproducibility
            
        Returns:
            X_train, X_test, y_train, y_test
        """
        X_train, X_test, y_train, y_test = train_test_split(
            features, target, test_size=test_size, random_state=random_state
        )
        
        return X_train, X_test, y_train, y_test
    
    def create_preprocessing_pipeline(
        self,
        numeric_features: List[str],
        categorical_features: List[str]
    ) -> ColumnTransformer:
        """
        Create a preprocessing pipeline for feature transformation.
        
        Args:
            numeric_features: List of numeric feature column names
            categorical_features: List of categorical feature column names
            
        Returns:
            ColumnTransformer preprocessing pipeline
        """
        numeric_transformer = Pipeline(steps=[
            ('imputer', SimpleImputer(strategy='median')),
            ('scaler', StandardScaler())
        ])
        
        categorical_transformer = Pipeline(steps=[
            ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
            ('onehot', OneHotEncoder(handle_unknown='ignore'))
        ])
        
        preprocessor = ColumnTransformer(transformers=[
            ('num', numeric_transformer, numeric_features),
            ('cat', categorical_transformer, categorical_features)
        ])
        
        return preprocessor
    
    def train_model(
        self,
        X_train: pd.DataFrame,
        y_train: pd.Series,
        model_type: str = 'random_forest',
        model_params: Optional[Dict] = None,
        cv: int = 5,
        feature_selection: bool = False,
        n_features_to_select: int = 10
    ) -> Tuple[Pipeline, Dict]:
        """
        Train a machine learning model.
        
        Args:
            X_train: Training features
            y_train: Training targets
            model_type: Type of model to train
            model_params: Parameters for the model
            cv: Number of cross-validation folds
            feature_selection: Whether to perform feature selection
            n_features_to_select: Number of features to select
            
        Returns:
            Trained model pipeline and cross-validation results
        """
        # Identify numeric and categorical features
        numeric_features = X_train.select_dtypes(include=['int64', 'float64']).columns.tolist()
        categorical_features = X_train.select_dtypes(include=['object', 'category']).columns.tolist()
        
        # Create preprocessing pipeline
        preprocessor = self.create_preprocessing_pipeline(numeric_features, categorical_features)
        
        # Create base model
        if model_type == 'random_forest':
            model = RandomForestClassifier(random_state=42)
            if model_params is None:
                model_params = {
                    'classifier__n_estimators': [100, 200],
                    'classifier__max_depth': [None, 10, 20],
                    'classifier__min_samples_split': [2, 5, 10]
                }
        elif model_type == 'gradient_boosting':
            model = GradientBoostingClassifier(random_state=42)
            if model_params is None:
                model_params = {
                    'classifier__n_estimators': [100, 200],
                    'classifier__learning_rate': [0.01, 0.1, 0.2],
                    'classifier__max_depth': [3, 5, 10]
                }
        elif model_type == 'svm':
            model = SVC(probability=True, random_state=42)
            if model_params is None:
                model_params = {
                    'classifier__C': [0.1, 1, 10],
                    'classifier__gamma': ['scale', 'auto'],
                    'classifier__kernel': ['rbf', 'linear']
                }
        elif model_type == 'logistic_regression':
            model = LogisticRegression(random_state=42)
            if model_params is None:
                model_params = {
                    'classifier__C': [0.1, 1, 10],
                    'classifier__penalty': ['l1', 'l2'],
                    'classifier__solver': ['liblinear', 'saga']
                }
        else:
            raise ValueError(f"Unsupported model type: {model_type}")
        
        # Create feature selection if requested
        if feature_selection:
            steps = [
                ('preprocessor', preprocessor),
                ('feature_selection', SelectKBest(f_classif, k=n_features_to_select)),
                ('classifier', model)
            ]
        else:
            steps = [
                ('preprocessor', preprocessor),
                ('classifier', model)
            ]
        
        # Create full pipeline
        pipeline = Pipeline(steps)
        
        # Perform grid search
        grid_search = GridSearchCV(
            pipeline,
            param_grid=model_params,
            cv=cv,
            scoring='roc_auc',
            n_jobs=-1
        )
        
        grid_search.fit(X_train, y_train)
        
        # Get best model
        best_model = grid_search.best_estimator_
        
        # Calculate cross-validation scores
        cv_scores = cross_val_score(
            best_model, X_train, y_train, cv=cv, scoring='roc_auc'
        )
        
        cv_results = {
            'mean_cv_score': cv_scores.mean(),
            'std_cv_score': cv_scores.std(),
            'best_params': grid_search.best_params_
        }
        
        # Store model and metrics
        model_key = f"{model_type}_{len(self.models)}"
        self.models[model_key] = best_model
        self.model_metrics[model_key] = {'cv_results': cv_results}
        
        return best_model, cv_results
    
    def evaluate_model(
        self,
        model: Pipeline,
        X_test: pd.DataFrame,
        y_test: pd.Series,
        model_key: Optional[str] = None
    ) -> Dict:
        """
        Evaluate a trained model on test data.
        
        Args:
            model: Trained model pipeline
            X_test: Test features
            y_test: Test targets
            model_key: Key to identify the model (optional)
            
        Returns:
            Dictionary of evaluation metrics
        """
        # Get predictions
        y_pred = model.predict(X_test)
        y_prob = model.predict_proba(X_test)[:, 1]
        
        # Calculate metrics
        metrics = {
            'accuracy': accuracy_score(y_test, y_pred),
            'precision': precision_score(y_test, y_pred),
            'recall': recall_score(y_test, y_pred),
            'f1': f1_score(y_test, y_pred),
            'roc_auc': roc_auc_score(y_test, y_prob)
        }
        
        # Store evaluation metrics
        if model_key is not None:
            if model_key in self.model_metrics:
                self.model_metrics[model_key]['test_metrics'] = metrics
            else:
                self.model_metrics[model_key] = {'test_metrics': metrics}
        
        return metrics
    
    def extract_feature_importance(
        self,
        model: Pipeline,
        X_train: pd.DataFrame,
        model_key: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Extract feature importances from a trained model.
        
        Args:
            model: Trained model pipeline
            X_train: Training features
            model_key: Key to identify the model (optional)
            
        Returns:
            DataFrame with feature importances
        """
        # Get feature names
        feature_names = X_train.columns.tolist()
        
        # Extract the classifier from the pipeline
        classifier = model.named_steps['classifier']
        
        # Extract feature importances (if available)
        if hasattr(classifier, 'feature_importances_'):
            importances = classifier.feature_importances_
        elif hasattr(classifier, 'coef_'):
            importances = classifier.coef_[0]
        else:
            self.logger.warning(f"Model does not have built-in feature importances attribute")
            return pd.DataFrame()
        
        # Create a DataFrame of feature importances
        importance_df = pd.DataFrame({
            'feature': feature_names,
            'importance': importances
        })
        
        # Sort by importance
        importance_df = importance_df.sort_values('importance', ascending=False)
        
        # Store feature importances
        if model_key is not None:
            self.feature_importances[model_key] = importance_df
        
        return importance_df
    
    def calculate_shap_values(
        self,
        model: Pipeline,
        X_sample: pd.DataFrame,
        model_key: Optional[str] = None,
        max_display: int = 20
    ) -> Dict:
        """
        Calculate SHAP values for model interpretability.
        
        Args:
            model: Trained model pipeline
            X_sample: Sample of features to calculate SHAP values for
            model_key: Key to identify the model (optional)
            max_display: Maximum number of features to display
            
        Returns:
            Dictionary with SHAP explainer and values
        """
        try:
            # Transform the input data
            preprocessor = model.named_steps['preprocessor']
            X_transformed = preprocessor.transform(X_sample)
            
            # Extract the classifier
            classifier = model.named_steps['classifier']
            
            # Create a SHAP explainer
            if hasattr(classifier, 'predict_proba'):
                explainer = shap.TreeExplainer(classifier) if hasattr(classifier, 'estimators_') else shap.KernelExplainer(
                    classifier.predict_proba, X_transformed
                )
            else:
                explainer = shap.KernelExplainer(classifier.predict, X_transformed)
            
            # Calculate SHAP values
            shap_values_result = explainer.shap_values(X_transformed)
            
            # Store SHAP values
            result = {
                'explainer': explainer,
                'shap_values': shap_values_result,
                'X_transformed': X_transformed
            }
            
            if model_key is not None:
                self.shap_values[model_key] = result
            
            return result
        
        except Exception as e:
            self.logger.error(f"Error calculating SHAP values: {e}")
            return {}
    
    def save_model(
        self,
        model_key: str,
        filepath: Optional[str] = None
    ) -> str:
        """
        Save a trained model to disk.
        
        Args:
            model_key: Key identifying the model to save
            filepath: Path to save the model to (optional)
            
        Returns:
            Path to the saved model
        """
        if model_key not in self.models:
            raise ValueError(f"Model key {model_key} not found")
        
        if filepath is None:
            # Use default path from config
            models_dir = self.config.get_output_path("models")
            os.makedirs(models_dir, exist_ok=True)
            filepath = os.path.join(models_dir, f"{model_key}.pkl")
        
        # Save the model
        with open(filepath, 'wb') as f:
            pickle.dump(self.models[model_key], f)
        
        # Also save metrics
        metrics_path = f"{os.path.splitext(filepath)[0]}_metrics.pkl"
        with open(metrics_path, 'wb') as f:
            pickle.dump(self.model_metrics.get(model_key, {}), f)
        
        return filepath
    
    def load_model(
        self,
        filepath: str,
        model_key: Optional[str] = None
    ) -> Pipeline:
        """
        Load a trained model from disk.
        
        Args:
            filepath: Path to the saved model
            model_key: Key to assign to the loaded model (optional)
            
        Returns:
            Loaded model
        """
        # Load the model
        with open(filepath, 'rb') as f:
            model = pickle.load(f)
        
        # Generate a key if not provided
        if model_key is None:
            model_key = f"loaded_model_{len(self.models)}"
        
        # Store the model
        self.models[model_key] = model
        
        # Try to load metrics
        metrics_path = f"{os.path.splitext(filepath)[0]}_metrics.pkl"
        if os.path.exists(metrics_path):
            with open(metrics_path, 'rb') as f:
                self.model_metrics[model_key] = pickle.load(f)
        
        return model


class OutcomePredictor:
    """Predict outcomes using trained machine learning models."""
    
    def __init__(self, config_obj=None):
        """
        Initialize the outcome predictor.
        
        Args:
            config_obj: Configuration object (default: use global config)
        """
        self.config = config_obj if config_obj is not None else config
        self.logger = logging.getLogger(__name__)
        self.feature_extractor = FeatureExtractor(config_obj)
        self.trainer = ModelTrainer(config_obj)
        self.models = {}
    
    def load_model(self, model_path: str, model_key: str) -> None:
        """
        Load a trained model for prediction.
        
        Args:
            model_path: Path to the saved model
            model_key: Key to identify the model
        """
        self.models[model_key] = self.trainer.load_model(model_path, model_key)
    
    def predict_treatment_response(
        self,
        variants_df: pd.DataFrame,
        clinical_data: Optional[pd.DataFrame] = None,
        model_key: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Predict treatment response based on genetic and clinical data.
        
        Args:
            variants_df: DataFrame with variant data
            clinical_data: DataFrame with clinical data (optional)
            model_key: Key of the model to use (optional)
            
        Returns:
            DataFrame with predictions
        """
        # Extract features
        variant_features = self.feature_extractor.extract_variant_features(variants_df)
        
        # Combine features
        features = self.feature_extractor.combine_features(variant_features, clinical_data=clinical_data)
        
        # Select model
        if model_key is None:
            # Use the first available model
            if not self.models:
                raise ValueError("No models available for prediction")
            model_key = next(iter(self.models))
        
        if model_key not in self.models:
            raise ValueError(f"Model key {model_key} not found")
        
        # Make predictions
        model = self.models[model_key]
        y_prob = model.predict_proba(features)
        y_pred = model.predict(features)
        
        # Create result DataFrame
        result = pd.DataFrame({
            'predicted_response': y_pred,
            'response_probability': y_prob[:, 1]
        })
        
        # Add sample identifier if available
        if 'sample_id' in variants_df.columns:
            result['sample_id'] = variants_df['sample_id']
        
        return result
    
    def predict_injury_risk(
        self,
        variants_df: pd.DataFrame,
        clinical_data: Optional[pd.DataFrame] = None,
        model_key: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Predict risk of radiation or drug-induced injury.
        
        Args:
            variants_df: DataFrame with variant data
            clinical_data: DataFrame with clinical data (optional)
            model_key: Key of the model to use (optional)
            
        Returns:
            DataFrame with predictions
        """
        # This is similar to predict_treatment_response but uses a different model
        # Extract features
        variant_features = self.feature_extractor.extract_variant_features(variants_df)
        
        # Combine features
        features = self.feature_extractor.combine_features(variant_features, clinical_data=clinical_data)
        
        # Select model
        if model_key is None:
            # Use the first available model
            if not self.models:
                raise ValueError("No models available for prediction")
            model_key = next(iter(self.models))
        
        if model_key not in self.models:
            raise ValueError(f"Model key {model_key} not found")
        
        # Make predictions
        model = self.models[model_key]
        y_prob = model.predict_proba(features)
        y_pred = model.predict(features)
        
        # Create result DataFrame
        result = pd.DataFrame({
            'predicted_injury': y_pred,
            'injury_probability': y_prob[:, 1]
        })
        
        # Add sample identifier if available
        if 'sample_id' in variants_df.columns:
            result['sample_id'] = variants_df['sample_id']
        
        return result


# Visualizations for model evaluation
def plot_roc_curve(model, X_test, y_test, title="ROC Curve"):
    """Plot ROC curve for model evaluation."""
    y_prob = model.predict_proba(X_test)[:, 1]
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    auc = roc_auc_score(y_test, y_prob)
    
    plt.figure(figsize=(10, 8))
    plt.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc='lower right')
    
    return plt.gcf()

def plot_feature_importance(importance_df, top_n=20, title="Feature Importance"):
    """Plot feature importance for model interpretability."""
    plt.figure(figsize=(12, 10))
    top_features = importance_df.head(top_n)
    sns.barplot(x='importance', y='feature', data=top_features)
    plt.title(title)
    plt.xlabel('Importance')
    plt.ylabel('Feature')
    plt.tight_layout()
    
    return plt.gcf()

def plot_confusion_matrix(y_true, y_pred, labels=None, title="Confusion Matrix"):
    """Plot confusion matrix for model evaluation."""
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(8, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
               xticklabels=labels, yticklabels=labels)
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title(title)
    
    return plt.gcf()

def plot_shap_summary(shap_values, X, max_display=20):
    """Plot SHAP summary for model interpretability."""
    plt.figure(figsize=(12, 10))
    shap.summary_plot(shap_values, X, max_display=max_display, show=False)
    plt.tight_layout()
    
    return plt.gcf()