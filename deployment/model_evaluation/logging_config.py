"""
Structured logging configuration for evaluator module.

Provides JSON-formatted logging for machine parsing alongside
human-readable console output.
"""

import logging
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Any


class JSONFormatter(logging.Formatter):
    """
    JSON-structured logging formatter for machine parsing.
    
    Outputs log records as JSON objects with:
    - timestamp: ISO 8601 format
    - level: Log level name
    - logger: Logger name
    - message: Log message
    - module: Source module
    - function: Source function
    - line: Source line number
    """
    
    def __init__(self, include_extra: bool = True):
        super().__init__()
        self.include_extra = include_extra
    
    def format(self, record: logging.LogRecord) -> str:
        log_obj = {
            'timestamp': datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z'),
            'level': record.levelname,
            'logger': record.name,
            'message': record.getMessage(),
            'module': record.module,
            'function': record.funcName,
            'line': record.lineno,
        }
        
        if self.include_extra and hasattr(record, 'extra_data'):
            log_obj['extra'] = record.extra_data
        
        if record.exc_info:
            log_obj['exception'] = self.formatException(record.exc_info)
        
        return json.dumps(log_obj, default=str)


class StructuredLogger:
    """
    Wrapper around logging.Logger with structured logging support.
    """
    
    def __init__(self, name: str, logger: Optional[logging.Logger] = None):
        self._logger = logger or logging.getLogger(name)
        self.name = name
    
    def _log(self, level: int, message: str, extra_data: Optional[dict] = None, exc_info: bool = False):
        if extra_data:
            extra = {'extra_data': extra_data}
            self._logger.log(level, message, extra=extra, exc_info=exc_info)
        else:
            self._logger.log(level, message, exc_info=exc_info)
    
    def debug(self, message: str, **kwargs):
        self._log(logging.DEBUG, message, kwargs or None)
    
    def info(self, message: str, **kwargs):
        self._log(logging.INFO, message, kwargs or None)
    
    def warning(self, message: str, **kwargs):
        self._log(logging.WARNING, message, kwargs or None)
    
    def error(self, message: str, **kwargs):
        self._log(logging.ERROR, message, kwargs or None, exc_info=kwargs.get('exc_info', False))
    
    def critical(self, message: str, **kwargs):
        self._log(logging.CRITICAL, message, kwargs or None, exc_info=kwargs.get('exc_info', False))


def setup_logging(
    log_file: Optional[Path] = None,
    console_level: int = logging.INFO,
    file_level: int = logging.DEBUG,
    log_format: str = 'json',
    log_dir: Optional[Path] = None,
) -> logging.Logger:
    """
    Configure structured logging for the evaluator.
    
    Args:
        log_file: Path to log file. If None, defaults to {log_dir}/evaluate.log
        console_level: Logging level for console output
        file_level: Logging level for file output
        log_format: 'json' for JSON, 'text' for human-readable
        log_dir: Directory for log files (used if log_file is None)
    
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    logger.handlers.clear()
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level)
    console_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    if log_file is None and log_dir is not None:
        log_file = Path(log_dir) / "evaluate.log"
    
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_file, mode='a')
        file_handler.setLevel(file_level)
        
        if log_format == 'json':
            file_handler.setFormatter(JSONFormatter())
        else:
            file_formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(funcName)s:%(lineno)d - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            file_handler.setFormatter(file_formatter)
        
        logger.addHandler(file_handler)
    
    logger.propagate = False
    
    return logger


def get_logger(name: str) -> StructuredLogger:
    """
    Get a structured logger instance.
    
    Args:
        name: Logger name (typically __name__)
    
    Returns:
        StructuredLogger instance
    """
    return StructuredLogger(name)


class LogContext:
    """Context manager for temporary log level changes."""
    
    def __init__(self, logger: logging.Logger, level: int):
        self.logger = logger
        self.level = level
        self.original_level = None
    
    def __enter__(self):
        self.original_level = self.logger.level
        self.logger.setLevel(self.level)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logger.setLevel(self.original_level)


def log_function_call(logger: logging.Logger):
    """
    Decorator to log function calls with arguments and results.
    
    Usage:
        @log_function_call(logger)
        def my_function(arg1, arg2):
            ...
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            logger.debug(
                f"Calling {func.__name__}",
                extra_data={'function': func.__name__, 'args': str(args), 'kwargs': str(kwargs)}
            )
            try:
                result = func(*args, **kwargs)
                logger.debug(
                    f"Completed {func.__name__}",
                    extra_data={'function': func.__name__, 'status': 'success'}
                )
                return result
            except Exception as e:
                logger.error(
                    f"Failed {func.__name__}: {str(e)}",
                    extra_data={'function': func.__name__, 'status': 'error', 'error': str(e)},
                    exc_info=True
                )
                raise
        return wrapper
    return decorator
